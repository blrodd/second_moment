classdef Origin

    properties
        databasename
        maindb
        orid
        origin
        event
        netmag
        mt
        egf_subset
        joinevent
        subset
        joinmt
        eqinfo
        stations
        arrivals
    end

    methods
        function EQ = Origin(db, orid)
            EQ.maindb = db;
            EQ.origin = dblookup_table(db, 'origin'); 
            
            EQ.event = dblookup_table(db, 'event');
            EQ.netmag = dblookup_table(db, 'netmag');
            EQ.mt = dblookup_table(db, 'mt');
            
            EQ.joinevent = dbjoin(EQ.origin, EQ.event);
            EQ.joinevent = dbjoin(EQ.joinevent, EQ.netmag);
        
            string = strcat('orid=="', num2str(orid),'"');
            table = dbsubset(EQ.joinevent, string);

            if dbnrecs(table) == 0
                elog_die(sprintf('Table subset for orid %s has no records', orid))
            end

            if dbnrecs(table) > 1
                % need to write something to only run if 1 result
                subset = dbsubset(table, 'netmag.auth!~/mt.SOCAL_MODEL/ && netmag.auth!="orbevproc"');
                [etime,elat,elon,edepth,mag,eauth,eorid,eevid]=dbgetv(subset, 'time', 'lat', 'lon', 'depth', 'magnitude', 'auth', 'orid', 'evid');
            end

            if dbnrecs(table) == 1
                [etime,elat,elon,edepth,mag,eauth,eorid,eevid]=dbgetv(table, 'time', 'lat', 'lon', 'depth', 'magnitude', 'auth', 'orid', 'evid');
            end

            % GET MOMENT TENSOR INFORMATION %
            subset = dbsubset(EQ.origin, string);
            mt_join = dbjoin(subset, EQ.mt);
            if dbnrecs(mt_join) > 0
                [estatus,emag,estrike1,estrike2,edip1,edip2]=dbgetv(mt_join, 'estatus', 'drmag', 'str1', 'str2', 'dip1', 'dip2');
                if (estatus == 'Quality: 0') | (estatus == 'Quality: 1')
                    mt_flag = 1;
                    elog_notify('MT Quality < 2: Do not use fault dimensions')
                else
                    mt_flag = 0;
                    elog_notify('MT Quality >= 2: Use fault dimensions')
                end
            else
              mt_flag = 1;
              elog_notify('MT solution does not exist')
            end

            if mt_flag == 1;
                strike1 = 307; dip1 = 83;
                strike2 = 216; dip2 = 82; % from moment tensor solution
                mag = mag;
                estatus = 'Quality: NA';
            else
                strike1 = estrike1; strike2 = estrike2;
                dip1 = edip1; dip2 = edip2;
                mag = emag;
            end

            EQ.eqinfo = struct('etime',etime, 'elat', elat, 'elon', elon, 'edepth', edepth, 'mag', mag, 'eauth', eauth, 'eorid', eorid, 'eevid', eevid, 'estatus',estatus, 'strike1', strike1, 'strike2', strike2, 'dip1', dip1, 'dip2', dip2);
        end

    
        function orids = getEGF(EQ, loc_margin, dep_margin, time_margin)
            string = strcat('evid !=', num2str(EQ.eqinfo.eevid), ' && magnitude < ', num2str(EQ.eqinfo.mag), ' && lon > ', num2str(EQ.eqinfo.elon-loc_margin), ...
                    ' && lon < ', num2str(EQ.eqinfo.elon+loc_margin), ' && lat > ', num2str(EQ.eqinfo.elat-loc_margin), ...
                    ' && lat < ', num2str(EQ.eqinfo.elat+loc_margin), ' && depth > ', num2str(EQ.eqinfo.edepth-dep_margin), ...
                    ' && depth < ', num2str(EQ.eqinfo.edepth+dep_margin), ' && time < ', num2str(EQ.eqinfo.etime+time_margin));
            egf_subset = dbsubset(EQ.joinevent, string);

            if dbnrecs(egf_subset) > 0
                [aorids,aevids,aprefors]=dbgetv(egf_subset, 'orid','evid', 'prefor');
                evids = unique(aevids);
                egf_orids = [];
          
                for i = length(evids);
                    evid = evids(i);
                    ind = find(aevids == evid);
                    orids = aorids(ind);
                    prefor = unique(aprefors(ind));
        
                    if find(orids == prefor) > 0
                        egforids = [ egf_orids ; prefor ];
                    else
                        egforids = [ egf_orids; orids(1) ];
                    end
                end
            orids = egforids;
            else
              elog_die(sprintf('No aftershock in database for orid %s', EQ.eqinfo.eorid))
            end
        end
        
        function EQ = get_stations(EQ, select, reject)
            MSyearday = epoch2str(EQ.eqinfo.etime, '%Y%j');
%            EGFyearday = epoch2str(EQ.EGFinfo.etime(index), '%Y%j')

            steps = {};
            steps{1} = 'dbopen site';
            
            
            steps{2} = sprintf('dbsubset ondate <= %s && (offdate >= %s || NULL)', MSyearday, MSyearday); 
            steps{3} = 'dbsort sta';

            if select
                steps{4} = sprintf('dbsubset sta=~/%s/', select);
            end

            if reject
                steps{5} = sprintf('dbsubset sta!~/%s/', reject);
            end
    
            
            
            site_table = dbprocess(EQ.maindb, steps);
            [stas, lats, lons, elevs] = dbgetv(site_table, 'sta', 'lat', 'lon', 'elev');
            EQ.stations = struct('sta', stas, 'lat', [], 'lon', [], 'elev', [], 'esaz', [], 'delta', [], 'distance', [], 'pdelay', [] ...
        ,'sdelay', [], 'pstime', []);

            for i=1:length(stas)
                site_table.record = i-1;

                esaz = dbeval(site_table, sprintf('azimuth(%s,%s,%s,%s)', ...
                            EQ.eqinfo.elat, EQ.eqinfo.elon, lats(i), lons(i))); 
                 
                delta = dbeval(site_table, sprintf('distance(%s,%s,%s,%s)', ...
                            EQ.eqinfo.elat, EQ.eqinfo.elon, lats(i), lons(i))); 
                realdistance = dbeval(site_table, sprintf('deg2km(%0.4f)', delta));

                pdelay = dbeval(site_table, sprintf('ptime(%0.4f,%0.4f)', delta, EQ.eqinfo.edepth));
                sdelay = dbeval(site_table, sprintf('stime(%0.4f,%0.4f)', delta, EQ.eqinfo.edepth));
                pstime = sdelay - pdelay;
 
                EQ.stations(i).lat = lats(i);
                EQ.stations(i).lon = lons(i);
                EQ.stations(i).elev = elevs(i);
                EQ.stations(i).esaz = esaz;
                EQ.stations(i).delta = delta;
                EQ.stations(i).distance = realdistance;                
                EQ.stations(i).pdelay = pdelay;
                EQ.stations(i).sdelay = sdelay;
                EQ.stations(i).pstime = pstime;
                
            end
            
        end

        function EQ = get_arrivals(EQ, select, reject)
            steps = {};
            steps{1} = 'dbopen assoc';
            steps{2} = 'dbjoin arrival';
            
            if select
                steps{3} = sprintf('dbsubset sta=~/%s/', select);
            end

            if reject
                steps{4} = sprintf('dbsubset sta!~/%s/', reject);
            end

            table = dbprocess(EQ.maindb, steps);
            % subset for MS orid
            string = sprintf('orid == %s', num2str(EQ.eqinfo.eorid));
            subset = dbsubset(table, string);
            subset = dbsubset(subset, 'auth!="anza"');
            
            % grab MS arrival info
            [stas, chans, iphases, phases, times, deltimes, snrs, esazs, deltas] = dbgetv(subset, 'sta', 'chan', 'iphase', 'phase', 'time', 'deltim', 'snr', 'esaz', 'delta');
            % store EQ arrival info 
            EQ.arrivals = struct('sta', stas, 'chan', chans, 'iphase', iphases ...
                   , 'phase', phases, 'time', num2cell(times), 'deltime', num2cell(deltimes) ...
                   , 'snr', num2cell(snrs), 'esaz', num2cell(esazs), 'delta', num2cell(deltas), 'pstime', []); 
            % add pstime from EQ.stations
            for i=1:length({EQ.arrivals.sta})
                sta = EQ.arrivals(i).sta;
                EQ = getPStime(EQ, sta)


            %for i=1:length({EQ.arrivals.sta})
            %    s = EQ.arrivals(i).sta;
            %   
            %    pid = find(strcmp({EQ.stations.sta}, s) == 1);
            %    EQ.arrivals(i).pstime = EQ.stations(pid).pstime;
            %end


        end % arrivals function

        function EQ=getPStime(EQ, sta)
            match = find(strcmp({EQ.arrivals.sta}, sta) == 1);
            if length(match) > 1
                pstime = abs(EQ.arrivals(match(1)).time - EQ.arrivals(match(2)).time)
                w = find(strcmp({EQ.arrivals.sta}, sta) == 1 && strcmp({EQ.arrivals.iphase}, 'P') == 1)
                EQ.arrivals(w).pstime = pstime
            else
                pid = find(strcmp({EQ.stations.sta}, sta) == 1);
                EQ.arrivals(w).pred_pstime = EQ.stations(pid).pstime

    end % methods

end % class
