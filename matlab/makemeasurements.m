function [t2,DONE,STF,GFsv,dhatsv,datasv,Tsv,T1sv,epsv,epldsv,tpldsv,t0,t1,PhaseSv]=makemeasurements(velEGFa,velMSa,npMS,npEGF,dtsva,stasm,compm,phasem, timems,timeegf, duration, niter, misfit_criteria, pickt2);
ns=size(velEGFa,1);

global mode_run
% Will save the measurements in MEASUREMENTS.MAT

% Initialize variables
DONE=zeros(ns,1); logical(DONE);
Npts=2048;
t2(1:ns)=99999; t1(1:ns)=0; t0=t1; dtsv(1:ns)=0;
datasv=zeros(ns,Npts);
dtsv=zeros(ns,1);
dhatsv=zeros(ns,Npts);
GFsv=zeros(ns,Npts);
STF=zeros(ns,Npts);
STF_sm=zeros(size(STF));
Tsv=zeros(ns,1); 
T1sv=zeros(ns,1);
epsv=zeros(ns,1);
epldsv=zeros(ns,Npts);
tpldsv=zeros(ns,Npts);
Mf(1:ns) = 0;
dofilt=logical(1);   
dtsv=dtsva;


% Automated mode.
if ~mode_run.interactive
    for i=1:ns
        if mode_run.verbose
            elog_notify(sprintf('%s_%s ASTF CALCULATION', stasm{i}, compm{i}))
        end
        % if EGF and MS samples > 100
        if(npEGF(i)>=100 && npMS(i)>=100)
            % grab MS data
            velMS=velMSa(i,1:npMS(i)); velMS=velMS';
  
            % grab EGF data
            velEGF=velEGFa(i,1:npEGF(i)); velEGF=velEGF';
 
            % time from start of waveform
            tms=[1:npMS(i)]*dtsva(i);
            tegf=[1:npEGF(i)]*dtsva(i);
            doi=1;
        else
            doi=0;
        end
        
        if doi == 1;
            %whos tms tegf velMS velEGF
            
            % If debug_plot mode, plot EGF and MS data together.
            if mode_run.debug_plot
                figure
    
                plot(tms-tms(1),velMS/max(velMS));
                hold on
                title([stasm{i},'-',compm{i}])
 
                plot(tegf-tegf(1),velEGF/max(velEGF),'r');
                legend('Mainshock','EGF'); xlabel('Time (s)');
                k = waitforbuttonpress;
                close
            end


         
            % Get samples since waveform start of arrival (8 seconds after start). 
            t = 8./dtsv(i);
        
            % Set MS inversion window.
            % If P phase, set window to end 3 seconds after arrival or to predicted P-S time if less than 3 seconds.
            % Prevents S-wave contamination. 
            if strcmp(phasem(i),'P') == 1 & duration
                len = floor(duration(i));
                 if len > 3
                   len = 3;
                 end
            else
                len = 3;
            end
 
            samps_before = 50;
            samps_after = len./dtsv(i);
            tt2b = t + samps_after;
            tt1b = t - samps_before;
    
            if mode_run.debug_plot
                figure
                plot(velMS);
                hold on
                vline(tt1b,'green')
                hold on
                vline(tt2b, 'green')
                hold on
                vline(t, 'red') 
                title(strcat('MS: ', stasm{i}, '-', compm{i}))
                k = waitforbuttonpress;
                close
            end

            % Do not allow data window to be too long, longer than Npts defined above.
            if(tt2b-tt1b>=Npts)
                data=velMS(tt1b:tt1b+Npts-1);
            else
                data=velMS(tt1b:tt2b); 
            end
            clf
            
            if mode_run.debug_plot
                figure
                plot(data);
                hold on
                vline(50, 'red') 
                title(strcat('MS Window: ', stasm{i}, '-', compm{i}))
                k = waitforbuttonpress;
                close
            end
            % add an interactive face if you want to change pick
            % could be useful to have option to add picks to stations without useable arrivals -- big rewrite
            % discuss with Juan, basically would grab data for every station present
            % for loop through that in this file and if use=='y', use that arrival
            % if use=='n' & in interactive mode, have option to pick arrivals
            % DONT WORRY ABOUT THIS FOR NOW
         
            % T1 is arrival - 3 seconds since inversion window
            T1=samps_before - 3;
            
            % taper inversion window
            np=length(data);
            zz=taper(np,.1);
            data=data.*zz; 
          
            % add trailing 0s after data until Npts
            data(np+1:Npts)=0;
         
            % new data matrix
            datasv(i,1:Npts)=data';
         
            % np = length of inversion window
            % T1 = MS arrival - 3 seconds since inversion window
            np=np-T1;
         
            % tt1b = arrival time since window beginning
            % tt2b = arrival time + MS inversion window length - (arrival - 3) - 1
            tt1b = t;
            tt2b=tt1b+np-1;
            
            if mode_run.debug_plot
                figure
                plot(velEGF);
                hold on
                vline(tt1b, 'green')
                vline(tt2b, 'green') 
                title(strcat('EGF: ', stasm{i}, '-', compm{i}))
                k = waitforbuttonpress;
                close
            end
         
            % EGF data is cut to EGF arrival time to tt2b 
            % length(EGF) is length(MS) - T1     
            GF=velEGF(tt1b:tt2b);
         
            if mode_run.debug_plot
                figure
                plot(GF);
                title(strcat('EGF Window: ', stasm{i}, '-', compm{i}))
                k = waitforbuttonpress;
            end
            % do not taper EGF
            zz=taper(np,.1);
            zz(1:np/2)=1;
            
            % add trailing 0s after data until Npts
            GF(np+1:Npts)=0;
         
            % new EGF data matrix
            GFsv(i,1:Npts)=GF';
            
            % get RSTF
            % inputs MS data, EGF data, T1, number of iterations
            [f,dhat,T,eps,tpld,epld]=pld(data,GF,T1,niter);
         
            % f - RSTF (apparent source time function)
            % dhat - RSTF * EGF (fit to the data seismogram
            % T - duration pick
            % eps - misfit of T
            % tpld - misfit
            % epld - duration tradeoff
         
            % finds 2nd moment of RSTF (t2) and mean centroid time (t1)
            [t2(i),t1(i),t0(i)]=findt2(f,pickt2);
            dt=dtsv(i);
            t2(i)=t2(i)*dt*dt;
            STF(i,1:Npts)=f';
            dhatsv(i,1:Npts)=dhat';
            Tsv(i)=(T-T1)*dt;
            T1sv(i)=T1*dt;
            epsv(i)=eps;
            npld=length(tpld);
            epldsv(i,1:npld)=epld(1:npld);
            tpldsv(i,1:npld)=tpld(1:npld);
            [junk,ind]=min(abs(tpld-T));
            PhaseSv(i) = phasem{i}; 
            if mode_run.verbose
                elog_notify(sprintf('   Misfit: %s ', epld(ind)))
                elog_notify(sprintf('   Apparent Duration: %s s', num2str(2*sqrt(t2(i)),2)))    
                elog_notify(sprintf('   ASTF Moment: %s', num2str(sum(STF(i,:)),5)))
            end
        
            if mode_run.debug_plot 
                % plot the 4 graphs
                figure
                dt=dtsva(i);
                subplot(2,2,1)
             
                % plot MS data
                plot([1:length(data)]*dt,data/max(data),'k'); hold on;
             
                % plot EGF data
                plot(dt*t1(i)+[1:length(GF)]*dt,GF/max(GF),'r');
                %xlim([0,length(data)*dt])
                xlim([0 5]); ylim([-1.1 1.1]);
                xlabel('Time (s)')
                legend('Data','EGF');
                title([stasm{i},' ',compm{i}]);  
                
                % plot ASTF
                subplot(2,2,2)
                plot([1:length(STF(i,:))]*dt,STF(i,:)); hold on;
                xlabel('Time (s)')
                title(['ASTF, moment:',num2str(sum(STF(i,:)),5)]);
                plot(t0(i)*dt,STF(round(t0(i))),'*')
                ylim([0 1.05*max(STF(i,:))])
                text(.1,.8*max(STF(i,:)),['ApprDur:',num2str(2*sqrt(t2(i)),2),' s'])
                xlim([0 2.1*t0(i)*dt]);
                %xlim([0,3*T1sv(i)])
             
                % plot misfit     
                subplot(2,2,3)
                plot(tpld*dt,epld); hold on;
                xlabel('Time (s)');
                ylabel('Misfit');
                xlim([0 3*t0(i)*dt])
                plot(tpld(ind)*dt,epld(ind),'*')
                ylim([0 1])
                
                % plot seismogram fit
                subplot(2,2,4)
                plot([1:length(data)]*dt,data,'k')
                hold on
                plot([1:length(dhat)]*dt,dhat,'r')
                legend('Data','EGF*STF')
                title(['Seismogram Fit']);
                %xlim([0,length(data)*dt])
                xlim([0 5]);
                xlabel('Time (s)')
                
                k = waitforbuttonpress;
                close
                end

            Mf(i) = epld(ind); 
            if epld(ind) > misfit_criteria
                elog_notify(sprintf('   DO NOT USE %s_%s: Misfit %.2f > %.2f', stasm{i}, compm{i}, epld(ind), misfit_criteria))
            else
                DONE(i) = 1;
            end
        end; %if
    end %for loop

% Flag ones with same station twice, and select the better result.
% Only use results of apparent duration in 1 std. 
    gids = find(DONE == 1);
    ads = 2*sqrt(t2(gids));
    inds = find(ads > (mean(ads) + 1 * std(ads)) | ads < (mean(ads) - 1 * std(ads)));
    DONE(gids(inds)) = 0;

    good = find(DONE == 1);
    [items, i] = unique({stasm{good}});
    for i=1:length(items)
        ids = find(ismember({stasm{good}}, items(i)));
        if length(ids) > 1
            [M, I] = max(Mf(good(ids))); 
            DONE(good(ids(I))) = 0;
        end 
    end
    good = find(DONE == 1);

    if mode_run.verbose
        elog_notify(sprintf('Usable stations: %s', strjoin(strcat({stasm{good}}, '_', {compm{good}}))))
    end

%                  %
% INTERACTIVE MODE %
%                  %


% for each waveform
else
 for i=1:ns
 
  % if EGF and MS samples > 100
  if(npEGF(i)>=100 && npMS(i)>=100)
   % grab MS data
   velMS=velMSa(i,1:npMS(i)); velMS=velMS';
  
   % grab EGF data
   velEGF=velEGFa(i,1:npEGF(i)); velEGF=velEGF';
 
   % time from start of waveform
   tms=[1:npMS(i)]*dtsva(i);
   tegf=[1:npEGF(i)]*dtsva(i);
   doi=1;
  else
   doi=0;
  end
 
  while(doi);
     %whos tms tegf velMS velEGF
  
     % initiate figure
     figure
 
     % plot(time, scaled MS)
     plot(tms-tms(1),velMS/max(velMS));
     hold on
     title([stasm{i},'-',compm{i}])
 
     % plot(time, scaled EGF)
     plot(tegf-tegf(1),velEGF/max(velEGF),'r');
     legend('Mainshock','EGF'); xlabel('Time (s)');
 
     % DONE is 0 to start
     if(DONE(i)==1)
       ipick(i)=menu('THIS ONE IS DONE, REDO?','Yes','No');
     else
       ipick(i)=menu('NOT DONE, Go Ahead w/Decon?','Yes','No');
     end
     close
 
     % if Yes was selected for Decon do the following
     if(ipick(i) == 1) 
     
      % plot MS data
      figure
      plot(velMS)
      title(stasm{i})
 
      % select window to zoom
      % select at least 50 samples before desired arrivals
      % select ~ 5 seconds after arrival or when amplitude returned to near-zero, do not include later arrivals
      disp('Pick range to Zoom in to')
      [x y] = ginput(1);
      tt1b=round(x);
      [x y]=ginput(1);
      tt2b=round(x);
 
      % plot MS data from start and end zoom selection
      plot(velMS(tt1b:tt2b))
      title([stasm{i},' ',compm{i}])
 
      % select range to invert
      % 20-50 samples before first arrival, when ~ 0 amplitude
      % 3-5 seconds after arrival or near zero amplitude
      disp('Pick range to invert')
      [x y] = ginput(1);
      [x2 y2]=ginput(1);
      tt2b=tt1b+round(x2);
      tt1b=tt1b+round(x);
      
      % do not allow data window to be too long, longer than Npts defined above
      if(tt2b-tt1b>=Npts)
       data=velMS(tt1b:tt1b+Npts-1);
      else
       data=velMS(tt1b:tt2b); 
      end
      clf
 
      % plot inversion range data
      plot(data)
 
      % pick arrival time
      disp('Pick P or S-wave arrival time');
      [x y] = ginput(1);
      
      % T1 is arrival - 3 seconds since inversion window
      T1=round(x)-3;
 
      % taper inversion window
      np=length(data);
      zz=taper(np,.1);
      data=data.*zz; 
  
      % add trailing 0s after data until Npts
      data(np+1:Npts)=0;
 
      % new data matrix
      datasv(i,1:Npts)=data';
      close
 
      % plot EGF data
      figure
      plot(velEGF)
      title(strcat('EGF: ', stasm{i}, '-', compm{i}))
 
      % pick range to zoom, window tightly around arrival
      disp('Pick range to Zoom')
      [x y] = ginput(1);
      tt1b=round(x);
      [x y]=ginput(1);
      tt2b=round(x);
 
      % plot zoom selection
      plot(velEGF(tt1b:tt2b))
      title([stasm{i},' ',compm{i}])
 
      % pick arrival time
      disp('Pick P or S wave arrival time')
      [x y] = ginput(1);
      close
      
      % tt1b = arrival time since window beginning
      tt1b=round(x)+tt1b; 
      
      %np comes from Mainshock
      %keyboard
      
      % np = length of inversion window
      % T1 = MS arrival - 3 seconds since inversion window
      np=np-T1;
 
      % tt2b = arrival time + MS inversion window length - (arrival - 3) - 1
      tt2b=tt1b+np-1;
 
      % EGF data is cut to EGF arrival time to tt2b 
      % length(EGF) is length(MS) - T1     
      GF=velEGF(tt1b:tt2b);
 
    if mode_run.debug_plot
        figure
        plot(GF);
        title(strcat('EGF Window: ', stasm{i}, '-', compm{i}))
        k = waitforbuttonpress;
    end

      % do not taper EGF
      zz=taper(np,.1);
      zz(1:np/2)=1;
      %GF=GF.*zz';   % TAPER GF    !!!!! DON"T TAPER GF BECAUSE IT STARTS at 
 %				    %ARRIVAL TIME?
     
      % add trailing 0s after data until Npts
      GF(np+1:Npts)=0;
 
      % new EGF data matrix
      GFsv(i,1:Npts)=GF';
      
      
      % get RSTF
      % inputs MS data, EGF data, T1, number of iterations
      [f,dhat,T,eps,tpld,epld]=pld(data,GF,T1,niter);
 
      % f - RSTF (apparent source time function)
      % dhat - RSTF * EGF (fit to the data seismogram
      % T - duration pick
      % eps - misfit of T
      % tpld - misfit
      % epld - duration tradeoff
 
      % finds 2nd moment of RSTF (t2) and mean centroid time (t1)
      [t2(i),t1(i),t0(i)]=findt2(f,pickt2);
      dt=dtsv(i);
      t2(i)=t2(i)*dt*dt;
      STF(i,1:Npts)=f';
      dhatsv(i,1:Npts)=dhat';
      Tsv(i)=(T-T1)*dt;
      T1sv(i)=T1*dt;
      epsv(i)=eps;
      npld=length(tpld);
      epldsv(i,1:npld)=epld(1:npld);
      tpldsv(i,1:npld)=tpld(1:npld);
      
      %keyboard
 
      % plot the 4 graphs
      figure
      dt=dtsva(i);
      subplot(2,2,1)
 
      % plot MS data
      plot([1:length(data)]*dt,data/max(data),'k'); hold on;
 
      % plot EGF data
      plot(dt*t1(i)+[1:length(GF)]*dt,GF/max(GF),'r');
      %xlim([0,length(data)*dt])
      xlim([0 5]); ylim([-1.1 1.1]);
      xlabel('Time (s)')
      legend('Data','EGF');
      title([stasm{i},' ',compm{i}]);  
      
      % plot ASTF
      subplot(2,2,2)
      plot([1:length(STF(i,:))]*dt,STF(i,:)); hold on;
      xlabel('Time (s)')
      title(['ASTF, moment:',num2str(sum(STF(i,:)),5)]);
      plot(t0(i)*dt,STF(round(t0(i))),'*')
      ylim([0 1.05*max(STF(i,:))])
      text(.1,.8*max(STF(i,:)),['ApprDur:',num2str(2*sqrt(t2(i)),2),' s'])
      xlim([0 2.1*t0(i)*dt]);
      %xlim([0,3*T1sv(i)])
 
      % plot misfit     
      subplot(2,2,3)
      plot(tpld*dt,epld); hold on;
      xlabel('Time (s)');
      ylabel('Misfit');
      xlim([0 3*t0(i)*dt])
      [junk,ind]=min(abs(tpld-T));
      plot(tpld(ind)*dt,epld(ind),'*')
      ylim([0 1])
      
      % plot seismogram fit
      subplot(2,2,4)
      plot([1:length(data)]*dt,data,'k')
      hold on
      plot([1:length(dhat)]*dt,dhat,'r')
      legend('Data','EGF*STF')
      title(['Seismogram Fit']);
      %xlim([0,length(data)*dt])
      xlim([0 5]);
      xlabel('Time (s)')
      
      % decide whether result is worth saving
      savei=menu('SAVE THIS RESULT','Yes P-wave','Yes S-wave', 'No REDO IT','No done');
      if(savei==1)
          DONE(i)=1;
          PhaseSv(i)='P';
          doi=0;
      elseif(savei==2);
          DONE(i)=1;
          PhaseSv(i)='S';
          doi=0;
      elseif(savei==3)
          DONE(i)=0;
          doi=1;
      elseif(savei==4)
          DONE(i)=0;
          doi=0;
      end
     elseif(ipick(i)==2)
      doi=0;
     end %ipick
     
  end; %while
  
  %%%% at the end of each station, save progress
  disp(['Done with ',stasm{i},' ',compm{i},' Saving to MEASUREMENTS.mat']);
  save('MEASUREMENTS.mat','DONE','Npts','t2','t1','t0','datasv','dtsv','dhatsv','GFsv',...
      'STF','STF_sm','Tsv','T1sv','epsv','epldsv','tpldsv');
 
 end  % stations
end
disp('Done Making Measurements');


