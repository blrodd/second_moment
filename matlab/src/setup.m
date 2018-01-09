%% runSJFsetup.m   
%% This is an example of how to do the setup if 
%% you have an antelope database and directories of
%% miniseed waveform files.
%% The goal is simply to populate the following list of variables,
%% you can replace this with any scheme you like
%% The idea is that you work of a desired list of staiton/components
%%
%% Output
%% compm  (Nsta,3) char variable of component names
%% depe   mainshock depth
%% dip1   dip of nodal plane #1
%% dip2   dip of nodal plane #2
%% dtsv   variable of individual component sample rates in seconds
%% late   mainshock latitude
%% lone   mainshock longitude
%% npEGF  number of points in each EGF waveform
%% npMS   number of points in each Mainshock waveform
%% slat   station latitudes
%% slon   station longitudes
%% stasm  cell array of station names
%% strike1 strike of nodal plane #1
%% strike2 strike of nodeal plane #2
%% velEGF  array of velocity seismograms for EGF
%% velMS   array of velocity seismograms for Mainshock

function [compm depe dip1 dip2 dtsv late lone npEGF npMS slat slon stasm strike1 strike2 velEGF velMS] = setup(verbose, debug, debug_plot, databasename, orid, select, reject, filter);

run('/opt/antelope/5.7/setup.m')

% Events and PLD parameters
dofilt=1;
fmin=0.25; fmax=15;  


%verbose = setflag('run_verbose')
%debug = setflag('run_debug')
%debug_plot = setflag('run_debug_plot')

%fprintf('verbose %f', verbose)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get origin info

% DO THIS IN  A CLASS STRUCTURE
% call arrivals class within it, to do del-times, etc.
% maybe call stations class within it as well, look at dbmoment

db = dbopen(databasename, 'r');
origin_table = dblookup_table(db, 'origin');
event_table = dblookup_table(db, 'event');
origin_event_table = dbjoin(origin_table, event_table);

string = strcat('orid=="', num2str(orid),'"');
origin_subset = dbsubset(origin_event_table, string);

[etime,elat,elon,edepth,eml,eauth,eorid,eevid]=dbgetv(origin_subset,'time','lat','lon','depth','ml','auth','orid','evid');  
%elat=str2num(elat); elon=str2num(elon); edepth=str2num(edepth);

mt_table = dblookup_table(db, 'mt');
join_table = dbjoin(origin_subset, mt_table);

if dbnrecs(join_table) > 0
    [estatus,emag,estrike1,estrike2,edip1,edip2]=dbgetv(join_table, 'estatus', 'drmag', 'str1', 'str2', 'dip1', 'dip2');
    
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

% paper indicated a mode you can run without mt solution, look into this.
% for example, leave strike/dip
% future - possible kill program if no fault info
if mt_flag == 1;
  strike1 = 307; dip1 = 83;
  strike2 = 216; dip2 = 82; % from moment tensor solution
  mag = eml;  
else
  strike1 = estrike1; strike2 = estrike2;
  dip1 = edip1; dip2 = edip2;
  mag = emag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get aftershock 
loc_margin = 0.002;
dep_margin = 10;
time_margin = 5*24*60*60; %5 day margin

string = strcat('evid !=', num2str(eevid), ' && ml < ', num2str(emag), ' && lon > ', num2str(elon-loc_margin), ...
        ' && lon < ', num2str(elon+loc_margin), ' && lat > ', num2str(elat-loc_margin), ...
        ' && lat < ', num2str(elat+loc_margin), ' && depth > ', num2str(edepth-dep_margin), ...
        ' && depth < ', num2str(edepth+dep_margin), ' && time < ', num2str(etime+time_margin));

egf_subset = dbsubset(origin_event_table, string);

if dbnrecs(egf_subset)>0
  [atimes,alats,alons,adepths,amls,aauths,aorids,aevids,aprefors]=dbgetv(egf_subset,'time','lat','lon','depth','ml','auth','orid','evid', 'prefor');  

  evids = unique(aevids);
  egf_orids = [];
  for i = length(evids);
    evid = evids(i);
    ind = find(aevids == evid);
    orids = aorids(ind);
    prefor = unique(aprefors(ind));
    if find(orids == prefor) > 0
      egf_orids = [ egf_orids ; prefor ];

    else
      egf_orids = [ egf_orids; orids(1) ];
    end
  end 

else
  elog_die(fprintf('No aftershock in database for orid %s', eorid))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Site info

%% change this to field structures
% need to connect this to arrivals table and only include stations with arrivals for the origin
site_table = dblookup_table(db, 'site');
string = strcat('sta=~/', select, '/');
site_table = dbsubset(site_table, string);

% grab info from this, sta, chan, slat, slon, selev, iphase

[sta, slat, slon, selev] = dbgetv(site_table, 'sta', 'lat', 'lon', 'elev');


str = 'Possible stations:';
for i = 1:length(sta)
  if i == 1
    str = strcat(str, '  ', sta{i});
  else 
    str = strcat(str, ',  ', sta{i});
  end
end

elog_notify(string)

% need to store event-station info
% do this in a class structure so you can apply functions to it. 
stations = [];
for i = 1:length(sta)
  stations.(sta{i}) = struct('lat', slat(i), 'lon', slon(i), 'elev', selev(i));
end

% look into Classes in matlab to do this instead of Structure, if possible
%station = struct('name', sta, 'lat', slat, 'lon', slon, 'elev', selev);
%disp(station)

assoc_table = dblookup_table(db, 'assoc')
arrival_table = dblookup_table(db, 'arrival')

site_table = dblookup_table(db, 'site');
string = strcat('sta=~/', select, '/');
site_table = dbsubset(site_table, string);

join_table = dbjoin(site_table, arrival_table);
join_table = dbjoin(join_table, assoc_table);
join_table = dbjoin(join_table, origin_subset);

[sta, chan, iphase, time, snr ] = dbgetv(site_table, 'sta', 'chan', 'iphase', 'time', 'snr');

%% do this in a class instead so you can apply functions to it such as calculating del-times, take-off angles, etc. 
arrivals = struct('sta', sta, 'chan', chan, 'iphase', iphase, 'time', time, 'snr', snr)

% use this to grab phase and channel code --> from there test on all channels, P on vertical, S on both horizontals.
% whichever produces best result of horizontal, use that. 

% need to set up recursiveness for find channel code similar to dbmoment -- Ask Juan about this

% once channel selected, then do waveform load for each station, channel 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figure
figure
plot(slon,slat,'k^','MarkerFaceColor','r','MarkerSize',10); hold on;
plot(elon,elat,'ko','MarkerFaceColor','b','MarkerSize',10);
axis('equal')

%keyboard
    
%indMS=find(eorid==MSorid); indMS=indMS(1);  % have more than one since not prefor
%MSyear = epoch2str(etimes, '%Y');
%MSyear=epoch2str(etimes(indMS),'%Y');

    
%db2=dbopen(inputdb,'r');
%db2=dblookup_table(db2,'origin');
%dum=strcat('orid=="',num2str(EGForid),'"');
%db2=dbsubset(db2,dum);
%[etimes,elats,elons,edepths,eauth,eorid,eevid]=dbgetv(db2,'time','lat','lon','depth','auth','orid','evid');   
%dbclose(db2);
%lategf=elats(1); lonegf=elons(1); depegf=edepths(1);
%indEGF=find(eorid==EGForid); indEGF=indEGF(1);  % have more than one since not prefor
%EGFyear=epoch2str(etimes(indEGF),'%Y');
%
%%Read in all waveforms
%%if from an accelerometer, then integrate to velocity
%npMS=zeros(ns,1);  npEGF=npMS;
%
%for i=1:ns
%    %i
%    zzf=strcat('*',char(stasm{i}),'*',compm(i,:),'*');    
%    if(compm(i,2)=='H')
%      dointe=0;
%    elseif(compm(i,2)=='N');
%      dointe=1;
%    end
%    
%    %MS datafile
%    dirname='MSdata';
%    %MS datafile
%    zz=['cd ',dirname];  %,num2str(MSevid)]
%    eval(zz);
%    f1=dir(zzf);  % find our file
%    file=f1.name;
%    X=rdmseed(file);
%    tms = cat(1,X.t);
%    dms = cat(1,X.d);
%    dtsv(i)=1./X(1).SampleRate;
%    sampr=X(1).SampleRate;
%    if(dofilt)
%       [B A] = butter(4,2*[fmin fmax]/sampr);
%       dms=dms-mean(dms);
%       dms=detrend(dms);
%       tp=taper(length(dms),.01);
%       dms=dms.*tp;
%       dms=filtfilt(B,A,dms);
%    end
%    npMS(i)=length(dms);
%
%    % convert to Velocity if needed
%    if(dointe)
%      velMS(i,1:npMS(i))=inte(dms-mean(dms),dtsv(i));
%    else
%      velMS(i,1:npMS(i))=dms-mean(dms);  
%    end
%    % Go back up and repeat for egf;
%    cd ../EGFdata   
%    f1=dir(zzf);  % find our file
%    
%    if(length(f1)==1)  % found right file    
%     file=f1.name;
%     X=rdmseed(file);
%     tegf = cat(1,X.t);
%     degf = cat(1,X.d);
%     if(dofilt)
%       degf=degf-mean(degf);
%       degf=detrend(degf);
%       tp=taper(length(degf),.05);
%       degf=degf.*tp;
%       degf=filtfilt(B,A,degf);
%     end
%     npEGF(i)=length(degf);
%     % convert to Velocity
%     if(dointe)
%      velEGF(i,1:npEGF(i))=inte(degf-mean(degf),dtsv(i));
%     else
%      velEGF(i,1:npEGF(i))=degf-mean(degf);
%     end
%    end %GF file exists
%    cd ..
%    
%end % loop over stations
%
%%keyboard
%save SJFsetup compm depe dip1 dip2 dtsv late lone npEGF npMS slat slon stasm strike1 strike2 velEGF velMS
%
%function notify(verbose, MSG)
%  if verbose == 1
%    elog_notify(MSG)
%  end
%
%function y = setflag(flag)
%  y = 0;
%  if exist(flag)
%    if strcmpi(eval(flag) , 'True')
%        y = 1;
%    end
%  end
%
