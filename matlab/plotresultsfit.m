%% script to plot the inversion result fit on a map
h = figure;
if ~mode_run.debug_plot
    set(h, 'visible', 'off')
end
 
subplot(1,2,1)
hh=scatter(slon(IJ),slat(IJ),100,2*sqrt(t2(IJ)),'filled'); hold on;
plot(slon(IJ),slat(IJ),'ko','MarkerSize',10);
colormap(jet); colorbar; caxis([0.15 0.50])
plot(lone,late,'w^','MarkerSize',15,'MarkerFaceColor','r');
axis([-117.0 -115.9 33.0 34.1])
xlabel('Longitude','FontSize',12);
ylabel('Latitude','FontSize',12);
set(gca,'FontSize',12);
title('Measurements \tau_c(s) in seconds ','FontSize',12);

subplot(1,2,2)
hh=scatter(slon(IJ),slat(IJ),100,2*sqrt(G*m2),'filled'); hold on;
plot(slon(IJ),slat(IJ),'ko','MarkerSize',10);
colormap(jet); colorbar; caxis([0.15 0.50])
plot(lone,late,'w^','MarkerSize',15,'MarkerFaceColor','r');
axis([-117.0 -115.9 33.0 34.1])
xlabel('Longitude','FontSize',12);
ylabel('Latitude','FontSize',12);
title('2nd Moments Fit, \tau_c(s) in seconds','FontSize',12);
set(gca,'FontSize',12);

if ~mode_run.no_figure
    saveas(h, sprintf('%s/MS%d_EGF%d_FITresult.png', image_dir, msorid, orid))
end

if mode_run.debug_plot
    k = waitforbuttonpress;
    close
end

