classdef Data
    properties
        time
        sta
        chan_code
        samprate
        ncalib
        segtype
        tr
        chan
        data
        dbview
    end %properties

    methods
        function WF = Data(db, sta, chan_code, time, tw, filter)
            % set starttime and endtime
            WF.time = time;
            starttime = WF.time;
            endtime = starttime + tw;
            
            WF.sta = sta;
            WF.chan_code = chan_code;
            %
            % GET WFDISC TABLE
            %
            steps = {};
            
            % open wfdisc
            steps{1} = 'dbopen wfdisc';

            % subset for station, chan_codenel, and time
            steps{2} = sprintf('dbsubset sta=~/%s/ && endtime > %s && time < %s && chan =~/%s/', ...
                                sta, starttime, endtime, chan_code); 
             % join to sensor, instrument
            steps{3} = 'dbjoin sensor';
            steps{4} = 'dbjoin instrument';
            steps{5} = 'dbsort sta chan';

            % set table view
            dbview = dbprocess(db, steps);

            % if no records return from function
            if dbnrecs(dbview) == 0
                elog_notify(sprintf('No traces after subset for sta=~/%s/ && chan=~/%s/', sta, chan_code))
                return
            end

            % if 3 records 
            if dbnrecs(dbview) == 3
                % get wfdisc info for later
                [samprate, ncalib,segtype, time] = dbgetv(dbview, 'wfdisc.samprate', ...
                                'instrument.ncalib', 'instrument.rsptype', 'wfdisc.time');
                %epoch2str(time(1), '%Y-%m-%d %H:%M:%S')
                %epoch2str(starttime, '%Y-%m-%d %H:%M:%S')
                
                WF.time = starttime - time(1);
                % only use one value for each, should be same for all 3 records
                WF.samprate = samprate(1); WF.ncalib = ncalib(1); WF.segtype = segtype(1); 

                %
                % LOAD DATA
                %
                
                % if data is not readable, return from function
                try
                    tr = trload_css(dbview, starttime, endtime);
                    %trsplice(tr)
                catch
                    elog_notify(sprintf('Could not read data for %s:%s', sta, chan_code))
                    return
                end
                
                % if no data, return from function
                if dbnrecs(tr) == 0
                    elog_notify(sprintf('No data after trload for %s:%s', sta, chan_code))
                    return
                end

                % if more than 3 records, return from function
                if dbnrecs(tr) > 3
                    elog_notify(sprintf('Too many traces after trload_cssgrp for %s:%s', sta, chan_code))
                    return
                end
                
                %          %
                % FIX DATA %
                %          %   

                % apply calibration
                % demean traces
                %trapply_calib(tr);
%                trfilter(tr, 'DEMEAN');
%
%                % apply integration
%                if strcmp(segtype, 'A')
%                    trfilter(tr, 'INT');
%                    segtype = 'V';
%
%                elseif strcmp(segtype, 'D')
%                    trfilter(tr, 'DIF');
%                    segtype = 'V';
%
%                elseif strcmp(segtype, 'V')
%                    segtype = 'V';
%
%                else
%                    elog_notify(sprintf('Unknown data type for %s:%s', sta, chan_code))
%                    return
%                end
%                
%                %  rotation -- discuss with Juan first???
%                % filter
%                try
%                    trfilter(tr, filter);
%                catch
%                    elog_notify(sprintf('Problems with filter %s for %s:%s', filter, sta, chan_code))
%                    return
%                end
%
%
                WF.tr = tr; 
                WF.dbview = dbview;
            end
        end % function

        function WF = grab_data(WF, chan, type)
            % subset tr for specific channel
            if strcmp(WF.segtype, 'A')
                dointe = 1;
            else
                dointe = 0;
            end
            
            dbview = dbsubset(WF.dbview, sprintf('chan =~/%s/', chan));
            [dir, file] = dbgetv(dbview, 'dir', 'dfile');
            filename = strcat(dir, '/', file);
            X=rdmseed(filename);
            tms = cat(1,X.t);
            dms = cat(1,X.d);
            dtsv=1./X(1).SampleRate;
            sampr=X(1).SampleRate;

            [B A] = butter(4, 2*[0.25 15]/sampr);
%            if (strcmp(chan, 'HNN') | strcmp(chan, 'HNZ') | strcmp(chan, 'HNE'))
%                [B A] = butter(4,2*[0.25 15]/sampr);
%            else 
%                [B A] = butter(4,2*[1 5]/sampr);
%            end
            dms=dms-mean(dms);

            dms=detrend(dms);
            
            if strcmp(type, 'ms')
                tp=taper(length(dms),.01);
            else
                tp=taper(length(dms),.05);
            end
            
            dms=dms.*tp;

            dms=filtfilt(B,A,dms);

            if(dointe)
              dms =inte(dms-mean(dms),dtsv);
            else
              dms =dms-mean(dms);
            end

            start = round(WF.time.*sampr);
             
            WF.data = dms(start:length(dms));

                % for longer mseed files, the tw will be changed
                % time: arrival - 8 seconds
                % original data before changes will be start: time - tw, end: time + 2*tw
                % there will be 2 extra time windows so filter and taper do not affect true window 
                % data = data[tw:2*tw]
                
            %elseif dbnrecs(tr) > 1
            %    elog_notify(sprintf('More than 1 trace for %s:%s', WF.sta, chan))
            %    return
            %elseif dbnrecs(tr) == 0
            %    elog_notify(sprintf('No traces for %s:%s', WF.sta, chan))
            %    return
            %end
        end % function 
            
    end % methods
end %class
