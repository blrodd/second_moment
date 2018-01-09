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


run('/opt/antelope/5.7/setup.m')
cd '/Users/rrodd/2ndMoment/2ndMomentsSuppMaterial_REV/SJFex'
inputdb='SJFexdb';
% Events and PLD parameters
niter=100;

% set parameters
% pull these from database and command line arguments in updated script
MSevid=1621;  emag=5.1;
EGFevid=1755; 
dointe=0;
dofilt=1;
fmin=0.25; fmax=15;  
strike1=307; dip1=83;
strike2=216; dip2=82; % from moment tensor solution  

% List of Components to do measurements on
stasm={'TRAN'; 'BALD'; 'IWR'; 'LVA2';  'SND';  'HSSP';   'BCCC'; 'CRY'; 'TRO'; 'SETM'; 'BZN'; 'BSAP'; 'RRSP'; 'MGD'; 'MGE'; 'SLR';'SAL'; 'DEV'; 'PSD'; 'MSC'; 'SNO'; 'PFO';};
% why does each station correspond with a different component???????
compm=['HNN';  'EHN';  'HNN'; 'HNE';   'HNE';   'HNE';   'HNE';  'HNN'; 'HNE'; 'HNN'; 'HHN'; 'HNE';   'HHN';  'HNE'; 'HNE'; 'HHN'; 'HHN'; 'HHZ'; 'HHZ'; 'HHZ'; 'HHZ'; 'HNZ';];

% Get info from antelope db
% ns = number of stations
ns=length(stasm); 
% for each station do the following:
for i=1:ns
  % open database
  db2=dbopen(inputdb,'r');
  % find site table
  db2=dblookup_table(db2,'site');  
  dum=strcat('sta=="',char(stasm{i}),'"');
  % subset for station
  db2=dbsubset(db2,dum);
  % get lat/lon for each station and save in vector
  [slat(i),slon(i)]=dbgetv(db2,'lat','lon');
  dbclose(db2);
end

% open and subset origin table, grab event info 
db2=dbopen(inputdb,'r');
db2=dblookup_table(db2,'origin');
dum=strcat('evid=="',num2str(MSevid),'"');
db2=dbsubset(db2,dum);
[etimes,elats,elons,edepths,eauth,eorid,eevid]=dbgetv(db2,'time','lat','lon','depth','auth','orid','evid');   
dbclose(db2);
late=elats(1); lone=elons(1); depe=edepths(1);

% plot stations and event info together
figure
plot(slon,slat,'k^','MarkerFaceColor','r','MarkerSize',10); hold on;
plot(lone,late,'ko','MarkerFaceColor','b','MarkerSize',10);
axis('equal')

%keyboard
% indMS = first evid in origin table    
indMS=find(eevid==MSevid); indMS=indMS(1);  % have more than one since not prefor
MSyear=epoch2str(etimes(indMS),'%Y');
    
% repeat for EGFevid
db2=dbopen(inputdb,'r');
db2=dblookup_table(db2,'origin');
dum=strcat('evid=="',num2str(EGFevid),'"');
db2=dbsubset(db2,dum);
[etimes,elats,elons,edepths,eauth,eorid,eevid]=dbgetv(db2,'time','lat','lon','depth','auth','orid','evid');   
dbclose(db2);
lategf=elats(1); lonegf=elons(1); depegf=edepths(1);
indEGF=find(eevid==EGFevid); indEGF=indEGF(1);  % have more than one since not prefor
EGFyear=epoch2str(etimes(indEGF),'%Y');

%Read in all waveforms
%if from an accelerometer, then integrate to velocity
% npMS = vector of length of each waveform)
npMS=zeros(ns,1);  npEGF=npMS;

% for each station
for i=1:ns
    %i
    zzf=strcat('*',char(stasm{i}),'*',compm(i,:),'*');    
    if(compm(i,2)=='H')
      dointe=0;
    elseif(compm(i,2)=='N');
      dointe=1;
    end
    
    %MS datafile
    dirname='MSdata';
    %MS datafile
    zz=['cd ',dirname];  %,num2str(MSevid)]
    eval(zz);
    f1=dir(zzf);  % find our file
    file=f1.name;
    X=rdmseed(file);
    tms = cat(1,X.t);
    dms = cat(1,X.d);
    % each waveform dt
    dtsv(i)=1./X(1).SampleRate;
    % each waveform sample rate
    sampr=X(1).SampleRate;
    % filter, demean, detrend, taper
    if(dofilt)
       [B A] = butter(4,2*[fmin fmax]/sampr);
       dms=dms-mean(dms);
       dms=detrend(dms);
       tp=taper(length(dms),.01);
       dms=dms.*tp;
       dms=filtfilt(B,A,dms);
    end
    
    % length of waveform saved for each in for loop
    npMS(i)=length(dms);

    % convert to Velocity if needed
    if(dointe)
      velMS(i,1:npMS(i))=inte(dms-mean(dms),dtsv(i));
    else
      % fill velMS for each waveform with 1: length of waveform minus the mean
      % seems to autofill remainder of vector of 0 till reaches max(npMS)
      % I do not see anywhere this is explicity coded. 
      velMS(i,1:npMS(i))=dms-mean(dms);  
    end
    % Go back up and repeat for egf;
    cd ../EGFdata   
    f1=dir(zzf);  % find our file
    
    if(length(f1)==1)  % found right file    
     file=f1.name;
     X=rdmseed(file);
     tegf = cat(1,X.t);
     degf = cat(1,X.d);
     if(dofilt)
       degf=degf-mean(degf);
       degf=detrend(degf);
       tp=taper(length(degf),.05);
       degf=degf.*tp;
       degf=filtfilt(B,A,degf);
     end
     npEGF(i)=length(degf);
     % convert to Velocity
     if(dointe)
      velEGF(i,1:npEGF(i))=inte(degf-mean(degf),dtsv(i));
     else
      velEGF(i,1:npEGF(i))=degf-mean(degf);
     end
    end %GF file exists
    cd ..
    
end % loop over stations

%keyboard
save SJFsetup compm depe dip1 dip2 dtsv late lone npEGF npMS slat slon stasm strike1 strike2 velEGF velMS
     
     
