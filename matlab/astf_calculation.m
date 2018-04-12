function [t2, done, stf,gfsv,dhatsv,datasv,tsv,t1sv,epsv,epldsv,tpldsv,t0,t1,phasesv] = astf_calculation(velMS, velEGF, dt, tms, tegf, sta, comp, phase, t, samps_before, samps_after, Npts, velMS_rot, velEGF_rot, niter, pickt2, misfit_criteria, update_arrival)
    global mode_run
    % set defaults for optional parameters
    if nargin < 17
        elog_error('Not enough input arguments in astf_calculation')
    elseif nargin == 17
        update_arrival = 0;
    elseif nargin > 18
        elog_error('Too many input arguments in astf_calculation')
    end  

    tt2b = t + samps_after;
    tt1b = t - samps_before;
    
    if update_arrival
        if strcmp(phase, 'P')
            xi = 0.5;
            nbins = [25 50 100 2/dt];
            s = 0;
        end

        if strcmp(phase, 'S')
            xi = 0.99;
            nbins = [10 15 20 50];
            s = 30;
        end
        [tt1b tt2b] = automated_arrival(xi, nbins, s, phase, 'MS', velMS_rot, tt1b, tt2b, samps_before, samps_after, dt);
        t = tt1b + samps_before;
    end

    if mode_run.debug_plot
        figure
        plot(velMS);
        hold on
        vline(tt1b,'green')
        hold on
        vline(tt2b, 'green')
        hold on
        vline(t, 'red')
        title(strcat('MS: ', sta, '-', comp))
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
        title(strcat('MS Window: ', sta, '-', comp))
        k = waitforbuttonpress;
        close
    end

    % T1 is arrival - 3 seconds since inversion window
    T1=samps_before - 3;

    % taper inversion window
    np=length(data);
    zz=taper(np,.1);
    data=data.*zz;

    % add trailing 0s after data until Npts
    data(np+1:Npts)=0;

    % new data matrix
    datasv(1:Npts)=data';

    % np = length of inversion window
    % T1 = MS arrival - 3 seconds since inversion window
    np=np-T1;

    % tt1b = arrival time since window beginning
    % tt2b = arrival time + MS inversion window length - (arrival - 3) - 1
    tt1b = t;
    tt2b=tt1b+np-1;
    
    if update_arrival
        if strcmp(phase, 'P')
            xi = 0.5;
            nbins = [25 50 100 2/dt];
            s = 0;
        end

        if strcmp(phase, 'S')
            xi = 0.99;
            nbins = [10 15 20 50];
            s = 30;
        end
        [tt1b tt2b] = automated_arrival(xi, nbins, s,  phase, 'EGF', velMS_rot, tt1b, tt2b, samps_before, samps_after, dt, np);
    end

    if mode_run.debug_plot
        figure
        plot(velEGF);
        hold on
        vline(tt1b, 'green')
        vline(tt2b, 'green')
        title(strcat('EGF: ', sta, '-', comp))
        k = waitforbuttonpress;
        close
    end

    % EGF data is cut to EGF arrival time to tt2b 
    % length(EGF) is length(MS) - T1     
    GF=velEGF(tt1b:tt2b);

    if mode_run.debug_plot
        figure
        plot(GF);
        title(strcat('EGF Window: ', sta, '-', comp))
        k = waitforbuttonpress;
    end
    % do not taper EGF
    zz=taper(np,.1);
    zz(1:np/2)=1;

    % add trailing 0s after data until Npts
    GF(np+1:Npts)=0;

    % new EGF data matrix
    gfsv(1:Npts)=GF';

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
    [t2,t1,t0]=findt2(f,pickt2);
    t2=t2*dt*dt;
    stf(1:Npts)=f';
    dhatsv(1:Npts)=dhat';
    tsv=(T-T1)*dt;
    t1sv=T1*dt;
    epsv=eps;
    npld=length(tpld);
    epldsv(1:npld)=epld(1:npld);
    tpldsv(1:npld)=tpld(1:npld);
    [junk,ind]=min(abs(tpld-T));
    phasesv = phase;
  
    if epld(ind) > misfit_criteria && ~update_arrival && mode_run.auto_arrival
        elog_notify(sprintf('Misfit %0.2f > Criteria %0.2f: Attempting to detect arrival to improve result', epld(ind), misfit_criteria)) 
        update_arrival = 1;
        [t2,done,stf,gfsv,dhatsv,datasv,tsv,t1sv,epsv,epldsv,tpldsv,t0,t1,phasesv] = astf_calculation(velMS, velEGF, dt, tms, tegf, sta, comp, phase, t, samps_before, samps_after, Npts, velMS_rot, velEGF_rot, niter, pickt2, misfit_criteria, update_arrival)    
    else
        if mode_run.verbose
            elog_notify(sprintf('   Misfit: %s ', epld(ind)))
            elog_notify(sprintf('   Apparent Duration: %s s', num2str(2*sqrt(t2),2)))
            elog_notify(sprintf('   ASTF Moment: %s', num2str(sum(stf(:)),5)))
        end

        if mode_run.debug_plot
            % plot the 4 graphs
            figure
            subplot(2,2,1)

            % plot MS data
            plot([1:length(data)]*dt,data/max(data),'k'); hold on;

            % plot EGF data
            plot(dt*t1+[1:length(GF)]*dt,GF/max(GF),'r');
            %xlim([0,length(data)*dt])
            xlim([0 5]); ylim([-1.1 1.1]);
            xlabel('Time (s)')
            legend('Data','EGF');
            title([sta,'-',comp]);

            % plot ASTF
            subplot(2,2,2)
            plot([1:length(stf(:))]*dt,stf(:)); hold on;
            xlabel('Time (s)')
            title(['ASTF, moment:',num2str(sum(stf(:)),5)]);
            plot(t0*dt,stf(round(t0)),'*')
            ylim([0 1.05*max(stf(:))])
            text(.1,.8*max(stf(:)),['ApprDur:',num2str(2*sqrt(t2),2),' s'])
            xlim([0 2.1*t0*dt]);
            %xlim([0,3*T1sv(i)])

            % plot misfit     
            subplot(2,2,3)
            plot(tpld*dt,epld); hold on;
            xlabel('Time (s)');
            ylabel('Misfit');
            xlim([0 3*t0*dt])
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

        if epld(ind) > misfit_criteria
            elog_notify(sprintf('   DO NOT USE %s_%s: Misfit %.2f > %.2f', sta, comp, epld(ind), misfit_criteria))
            done = 0;
        else
            done = 1;
        end
    end

