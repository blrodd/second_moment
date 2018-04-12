% runSJFex.m For given orid and set parameters, calculate the second moment and associated errors.
% Inputs:
% From second_moment.xpy
%   db:           database name given by first command-line argument
%   orid:         orid for main-shock event given by second command-line argument
%   select:       regex expression of stations to use
%   reject:       if not None, regex expression of stations to reject
%   filter:       butterworth bandpass filter for waveforms
%   tw:           time window in seconds for windowing waveforms around arrivals
%   loc_margin:   degrees that egf lat/lon can vary from mainshock lat/lon
%   dep_margin:   km that egf depth can vary from mainshock depth
%   time_margin:  time (s) that egf time can vary from mainshock time
% Outputs:
% Usage:
% Exceptions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiate Antelope-matlab interface.
run('/opt/antelope/5.8pre/setup.m')

% Add matlab files to MATLABPATH.
addpath(genpath('/Users/rrodd/bin/second_moment/matlab'));

% Modes to run data in (e.g. verbose, interactively, etc.)
global mode_run;

mode_run = struct('verbose', verbose, 'debug', debug, 'debug_plot', debug_plot, 'interactive', interactive, 'no_figure', no_figure, 'auto_arrival', auto_arrival);
fields = fieldnames(mode_run);
for i=1:numel(fields)
    if strcmp(mode_run.(fields{i}), 'False') == 1
        mode_run.(fields{i}) = false;
    else
        mode_run.(fields{i}) = true;
    end
end
% Add verbosity mode throughout code

% Get mainshock origin, station, arrival information.
db = dbopen(db, 'r');

% Initiate Origin class.
MS = Origin(db, orid);
MS = get_stations(MS, select, reject);
MS = get_arrivals(MS, select, reject);

% Set mainshock late, lone, depe for output.
late = MS.eqinfo.elat; lone = MS.eqinfo.elon; depe = MS.eqinfo.edepth;

% Set strikes and dips for output.
strike1 = MS.eqinfo.strike1; strike2 = MS.eqinfo.strike2;
dip1 =MS.eqinfo.dip1; dip2 = MS.eqinfo.dip2;

% Get all possible egf orids given the conditions.
msorid = orid;

% Determine how to input egf
if egf==-99
    egforids = getEGF(MS, loc_margin, dep_margin, time_margin);
else
    egforids = [egf];
end

% For each EGF possibility, run program.
EGFlen = length(egforids);
for index=1:EGFlen
    orid = egforids(index);
    filename = sprintf('%s/%d_setup_%d.mat', temp_dir, msorid, orid);

    % Add a data file flag if want to use specific data file 
    if LOADDATAFILE
        if mode_run.verbose
            elog_notify(sprintf('Load  waveforms from %s', filename))
        end

        try
            load filename
        catch
            elog_notify(sprintf('Filename %s does not exist.', filename))
            elog_die('Check to see if file exists or in correct location.\nRun with -e flag to specify egf orid')
        end            
    else
        if mode_run.verbose
            elog_notify(sprintf('Extract waveforms for MSorid %d / EGForid %d', msorid, orid))
        end
        tw = str2num(tw);
        [compm dtsv npEGF npMS slat slon stasm phasem velEGF velMS velEGF_rot velMS_rot timems timeegf duration] ...
                    = data_setup(db, orid, MS, select, reject, filter, tw);
        save filename compm depe dip1 dip2 dtsv late lone npEGF npMS slat slon ...
                    stasm phasem strike1 strike2 velEGF velMS timems timeegf duration 
    end

% Parameters for the PLD measurements of ASTF duration.
    measurements = sprintf('%s/%d_measurements_%d.mat', temp_dir, msorid, orid);
    if DOMEAS == 0
        if mode_run.verbose
            elog_notify(sprintf('Loading ASTF measurements from %s', measurements))
        end
        load measurements
    else
    % % call routine to make measurements
        if mode_run.verbose
            elog_notify(sprintf('Calculating ASTF measurements for MSorid %d / EGForid %d', msorid, orid))
        end
        [t2,DONE,STF,GFsv,dhatsv,datasv,T,T1sv,epsv,epld,tpld,t0,t1,PhaseSv] ...
                        = makemeasurements(velEGF,velMS,velEGF_rot, velMS_rot, npMS,npEGF,dtsv,stasm, ...
                        compm, phasem, timems,timeegf, duration, NITER, misfit_criteria, PICKt2);
        save measurements t2 DONE STF GFsv dhatsv datasv T T1sv epsv epld tpld t0 t1 PhaseSv
    end

    IJ=find(DONE==1);  %index of what stations had succesful measurements
%%keyboard
%
    %make map of results
    makeresultmap
    
    %plot the individual station summary plots
    makestationplots
end
%

%%keyboard
if(DOINVERSION)    
    % Define data vector for inversion; t2 values and phase=P or S
    IJ=find(DONE==1);
    clear d phas mlats mlons melevs
    for i=1:length(IJ)
      ii=IJ(i); 
      mlats(i)=slat(ii); mlons(i)=slon(ii); melevs(i)=0;  % FIX STRUCTURE selev(ii);
      d(i)=t2(ii);
      phas(i)=PhaseSv(ii);
    end

    % Define velocity model for partial calculation
    load Zigone2015v3;
    Vs=Vp/1.73;
    
    % Define fault plane
    strike=strike1; dip=dip1;
    
    % Get partials
    [G]=getpartials_2d_generic(mlats,mlons,melevs,late,lone,depe,Vp,Vs,topl,phas,strike,dip);
        
    % Finally do the inversion    
    [m2,L_c,W_c,vx,vy,tauc]=seconds_2d_v2(G,d');
    m2=m2(1:6);   % Ditch the dummy variable in the decision vector

    % Calculate some random things
    X=[m2(4), m2(5); m2(5), m2(6);];
    ssqr=sum((G*m2-d').^2)./sum(d.^2);  %useful definitions of misfit
    tauc_resids=sum(abs(2*sqrt(G*m2)-2*sqrt(d')));
    tauc_resids_norm=tauc_resids/sum(2*sqrt(d'));

    % Characteristic quantities
    [U,S,V]=svd(X);
    L_c=2*sqrt(max(max(S)));
    W_c=2*sqrt(S(2,2));
    v0=m2(2:3)/m2(1);
    mv0=sqrt(sum(v0.^2));
    tauc=2*sqrt(m2(1));
    L0=tauc*mv0;
    ratio=L0/L_c;
    ssqr3=sum((G*m2-d').^2)/sum((d'-mean(d)').^2);
    [m2', tauc,L_c,W_c,mv0,L0/L_c, ssqr3];

    plotresultsfit

    elog_notify(sprintf('tau_c %s (s), L_c %s (km), W_c %s (km), mv0 %s (km/s), L0/Lc %s', num2str(tauc,3), num2str(L_c,3), num2str(W_c,3), num2str(mv0,3), num2str(L0/L_c,3)));
end
%
%disp('Done with Inversion')
%disp('Type "return" to run Jackknife calculation')
%%keyboard
% if(DOJACKKNIFE)
%   clear mests L_cJK W_cJK vxJK vyJK taucJK
%   %if you only want the leave-one-out jacknife, this would suffice
%   %jackstat=jackknife(@seconds_2d_v2,G2,dd);     
%   %we want to delete all arrivals in a particular azimuth bin
%   %following McGuire, Zhao, Jordan, 2001 GJI.  
%   % get azimuths
%   clear dist az delta
%   ns=length(d);
%   for ii=1:ns
%     [dist(ii),az(ii),delta(ii)] = distaz(late,lone,mlats(ii),mlons(ii)); %lta,lna,ltb,lnb);
%   end
%   nband=360/azband;
%   for ii=1:nband
%      az1=(ii-1)*azband; az2=ii*azband;
%      IJ=find(az>=az1 & az<=az2);
%      Iuse=setdiff([1:ns],IJ);
%      G2ii=G(Iuse,:);   d2ii=d(Iuse);
%      [m2dum,L_cJK(ii),W_cJK(ii),vxJK(ii),vyJK(ii),taucJK(ii)]=seconds_2d_v2(G2ii,d2ii');
%      mests(ii,1:6)=m2dum(1:6);
%   end       
%   % Get Cij
%   clear xcm1bar varjack
%   xcm1bar(1:6)=0;
%   for i=1:nband
%     xcm1bar(1:6)=xcm1bar(1:6)+mests(ii,1:6);
%   end
%   xcm1bar=xcm1bar/nband;
%   varjack=zeros(6,6);
%   for ii=1:nband
%     clear dum;
%     dum=mests(ii,1:6) - xcm1bar;
%     varjack=varjack+dum'*dum;
%   end
%   varjack=varjack/(nband-1);
%   %calculates errors on the derived quantities tauc, Lc Wc v0, mv0, L0/Lc
%   [sigmatc,sigmaLc,sigmaWc,sigratio,sigv01,sigv02,sigmamV0]=geterrors(m2,varjack);     
%       
%   b=sprintf('%4.2f+-%3.3f  %4.2f+-%3.2f  %4.2f+-%3.2f  [%4.2f, %4.2f]+-[%3.2f, %3.2f]   %4.2f+-%3.2f  %4.2f+-%3.2f' ,...
%   tauc,sigmatc,L_c,sigmaLc,W_c,sigmaWc,v0(1),v0(2),sigv01,sigv02,mv0,sigmamV0,ratio,sigratio);
%   disp(['tauc          Lc            Wc              V0                     mV0           L0/LC'])
%   disp(b)
%        
%   %a=sprintf('%i %s %7.5f %8.5f %5.2f %2.1f %i %i %4.2f %5.3f %5.3f %5.3f %5.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f', ...
%   %MSevid,epoch2str(etimes(indMS),'%Y/%m/%d %H:%M:%S'), ... 
%   %late,lone,depe,emag, strike, dip, tauc, L_c, W_c, mv0, L0/L_c, ...
%   %m2(1), m2(2), m2(3), m2(4), m2(5), m2(6) ...
%   %)
% end   %End of Jackknife section
%
%disp('type return for bootstrap')
%%keyboard
%
%if(DOBOOTSTRAP)
%    
%   [mv0u,mv0l,bound2u,bound2l,Lcu,Lcl,taucu,taucl,minvrsv,minvru,minvrl]=bootstrap2nds(G,d,bconf,NB);
%   b=sprintf('                            [%4.2f,%4.2f],  [%4.2f,%4.2f],  [%4.2f,%4.2f]',taucl,taucu,Lcl,Lcu,mv0l,mv0u);
%   disp([num2str(bconf),'% confidence limits for tauc (s)       Lc (km)      mV0 (km/s)   '])
%   disp(b)
%
%end
% 
% 
