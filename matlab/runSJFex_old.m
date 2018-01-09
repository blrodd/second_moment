% change this to your 'src' directory
addpath('/Users/rrodd/bin/second_moment/src');

%cd '/Users/rrodd/bin/second_moment/SJFex'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flags and control parameters

LOADDATAFILE=0; %1 to just load the seismograms from an example .mat file
                %0 to read in miniseed files
PLOTMEAS=0;     %1 to make plots of the resulting STFs and data fits
DOMEAS=0;       %1 to go through the EGF deconvolutions station by station;
                %0 just loads example measurements .mat file
PICKt2=0;       %0 calculate mu^(0,2) (s) by integrating over whole ASTF
                %1 if you want to have add a second round of interactive
                %picking of the integration limits.
DOINVERSION=1;  %1 to run the convex optimization inversion


DOJACKKNIFE=1;  %1 to do jackknife error calculations.
azband=20;      %(degrees) size of azimuth bins to delete for jackknife error calcs.
DOBOOTSTRAP=1;  %1 to do bootstrap error calculations.
NB=1000;        % number of bootstrap iterations to run
bconf=0.95;     %confidence level to determine bootsrap bounds for
niter=100;      %number of iterations in the PLD deconvolution;  >50 is usually needed.
%strike,dip     % careful, these are set below after the earthquake data
%               % has been read in; they specify the plane for the 2nd moments
%               % inversion. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do your data input here;   
%%% Goal is to have the following arrays
%%% stasm   list of stations to try measurements on
%%% compm   associated list of individual components to try
%%% velMS   mainshock velocity seismograms, same order as stasm compm
%%% velEGF  EGF velocity seismograms, same order as stasm compm
%%% slats, slons   station latitude+longitudes
%%% eq
if(LOADDATAFILE)
  load SJFsetup.mat
else
  runSJFsetup
end

%%%%%%%%%%%%%Parameters for the PLD measurements of ASTF duration
if(DOMEAS)
 % call routine to make measurements
 [t2,DONE,STF,GFsv,dhatsv,datasv,T,T1sv,epsv,epld,tpld,t0,t1,PhaseSv]=makemeasurements(velEGF,velMS,niter,npMS,npEGF,dtsv,stasm,compm,PICKt2);
 IJ=find(DONE==1);  %index of what stations had succesful measurements
else
 load SJF_EXAMP_MEASUREMENTS
 IJ=find(DONE==1);
end
%keyboard

if(PLOTMEAS)
 %make map of results
 makeresultmap
 %plot the individual station summary plots
 makestationplots
end

clear zz;
for i=1:length(IJ)
    zz{i}=[char(stasm{IJ(i)}),' ',PhaseSv(IJ(i)),' tauc is ',num2str(2*sqrt(t2(IJ(i))),3),' and moment is=',num2str(sum(STF(IJ(i),:)),5)];
end
char(zz)
disp('DONE WITH MEASUREMENTS')
disp('Type "return" to run the inversion')
%keyboard
if(DOINVERSION)    
%%%%% Here is the Vanilla Inversion, just Volume>=0 constraint
%%
%% DEFINE DATA VECTOR FOR INVERSION; t2 values and phas=P or S
   IJ=find(DONE==1);
   clear d phas mlats mlons melevs
   for i=1:length(IJ)
     ii=IJ(i); 
     mlats(i)=slat(ii); mlons(i)=slon(ii); melevs(i)=0;  % FIX STRUCTURE selev(ii);
     d(i)=t2(ii);
     phas(i)=PhaseSv(ii);
   end

%% DEFINE VELOCITY MODEL FOR PARTIAL CALCULATION
    %Get Dimitri Zigone's BSSA 2015 model 
    %downsamled a lot for raytracer
    load Zigone2015v3;
    Vs=Vp/1.73;
    
%% define fault plane
    strike=strike1; dip=dip1;
    %strike=strike2; dip=dip2;   % pick a fault plane
    
% get partials
    [G]=getpartials_2d_generic(mlats,mlons,melevs,late,lone,depe,Vp,Vs,topl,phas,strike,dip);
    
% Finally do the inversion    
  [m2,L_c,W_c,vx,vy,tauc]=seconds_2d_v2(G,d');
  % m2=seconds_2d_v2(G,d2)
  m2=m2(1:6);   % Ditch the dummy variable in the decision vector
  % calculate some random things
  X=[m2(4), m2(5); m2(5), m2(6);];
  %X3=[m2(1), m2(2), m2(3); m2(2), m2(4), m2(5); m2(3), m2(5), m2(6)];
  ssqr=sum((G*m2-d').^2)./sum(d.^2);  %useful definitions of misfit
  tauc_resids=sum(abs(2*sqrt(G*m2)-2*sqrt(d')));
  tauc_resids_norm=tauc_resids/sum(2*sqrt(d'));
  % characteristic quantities
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
  %zz=[slat(K), slon(K), G*m2]
  %save LSPREDICTEDDATA -ASCII zz
  % MAKE STRING FOR RESULTS TABLE
  % EVID    Date     lat lon depth  Ml   strike dip tauc    Lc  Wc   mv0  DIR  tt xt yt xx xy yy
  % a=sprintf('%i %s %7.5f %8.5f %5.2f %2.1f %i %i %4.2f %5.3f %5.3f %5.3f %5.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f', ...
  % MSevid,epoch2str(etimes(indMS),'%Y/%m/%d %H:%M:%S'), ... 
  % elats(indMS),elons(indMS),edepths(indMS),emag, ...
  % strike, dip, tauc, L_c, W_c, mv0, L0/L_c, ...
  % m2(1), m2(2), m2(3), m2(4), m2(5), m2(6) ...
  % )
  plotresultsfit
  disp('tau_c (s), L_c (km), W_c (km), mv0 (km/s), L0/Lc are');
  disp([num2str(tauc,3),'    ',num2str(L_c,3),'    ',num2str(W_c,3),'    ',num2str(mv0,3),'    ',num2str(L0/L_c,3)]);
end

disp('Done with Inversion')
disp('Type "return" to run Jackknife calculation')
%keyboard
 if(DOJACKKNIFE)
   clear mests L_cJK W_cJK vxJK vyJK taucJK
   %if you only want the leave-one-out jacknife, this would suffice
   %jackstat=jackknife(@seconds_2d_v2,G2,dd);     
   %we want to delete all arrivals in a particular azimuth bin
   %following McGuire, Zhao, Jordan, 2001 GJI.  
   % get azimuths
   clear dist az delta
   ns=length(d);
   for ii=1:ns
     [dist(ii),az(ii),delta(ii)] = distaz(late,lone,mlats(ii),mlons(ii)); %lta,lna,ltb,lnb);
   end
   nband=360/azband;
   for ii=1:nband
      az1=(ii-1)*azband; az2=ii*azband;
      IJ=find(az>=az1 & az<=az2);
      Iuse=setdiff([1:ns],IJ);
      G2ii=G(Iuse,:);   d2ii=d(Iuse);
      [m2dum,L_cJK(ii),W_cJK(ii),vxJK(ii),vyJK(ii),taucJK(ii)]=seconds_2d_v2(G2ii,d2ii');
      mests(ii,1:6)=m2dum(1:6);
   end       
   % Get Cij
   clear xcm1bar varjack
   xcm1bar(1:6)=0;
   for i=1:nband
     xcm1bar(1:6)=xcm1bar(1:6)+mests(ii,1:6);
   end
   xcm1bar=xcm1bar/nband;
   varjack=zeros(6,6);
   for ii=1:nband
     clear dum;
     dum=mests(ii,1:6) - xcm1bar;
     varjack=varjack+dum'*dum;
   end
   varjack=varjack/(nband-1);
   %calculates errors on the derived quantities tauc, Lc Wc v0, mv0, L0/Lc
   [sigmatc,sigmaLc,sigmaWc,sigratio,sigv01,sigv02,sigmamV0]=geterrors(m2,varjack);     
       
   b=sprintf('%4.2f+-%3.3f  %4.2f+-%3.2f  %4.2f+-%3.2f  [%4.2f, %4.2f]+-[%3.2f, %3.2f]   %4.2f+-%3.2f  %4.2f+-%3.2f' ,...
   tauc,sigmatc,L_c,sigmaLc,W_c,sigmaWc,v0(1),v0(2),sigv01,sigv02,mv0,sigmamV0,ratio,sigratio);
   disp(['tauc          Lc            Wc              V0                     mV0           L0/LC'])
   disp(b)
        
   %a=sprintf('%i %s %7.5f %8.5f %5.2f %2.1f %i %i %4.2f %5.3f %5.3f %5.3f %5.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f', ...
   %MSevid,epoch2str(etimes(indMS),'%Y/%m/%d %H:%M:%S'), ... 
   %late,lone,depe,emag, strike, dip, tauc, L_c, W_c, mv0, L0/L_c, ...
   %m2(1), m2(2), m2(3), m2(4), m2(5), m2(6) ...
   %)
 end   %End of Jackknife section

disp('type return for bootstrap')
%keyboard

if(DOBOOTSTRAP)
    
   [mv0u,mv0l,bound2u,bound2l,Lcu,Lcl,taucu,taucl,minvrsv,minvru,minvrl]=bootstrap2nds(G,d,bconf,NB);
   b=sprintf('                            [%4.2f,%4.2f],  [%4.2f,%4.2f],  [%4.2f,%4.2f]',taucl,taucu,Lcl,Lcu,mv0l,mv0u);
   disp([num2str(bconf),'% confidence limits for tauc (s)       Lc (km)      mV0 (km/s)   '])
   disp(b)

end
 
 
