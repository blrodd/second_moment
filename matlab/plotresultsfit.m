%% script to plot the inversion result fit on a map
if mode_run.debug_plot
    figure('pos', [10 10 900 500]);
end

subplot(1,3,1)
hh=scatter(slon(IJ),slat(IJ),100,2*sqrt(t2(IJ)),'filled'); hold on;
plot(slon(IJ),slat(IJ),'ko','MarkerSize',10);
colormap(jet); colorbar; caxis([0.15 0.50])
plot(lone,late,'w^','MarkerSize',15,'MarkerFaceColor','r');
axis([-117.0 -115.9 33.0 34.1])
xlabel('Longitude','FontSize',12);
ylabel('Latitude','FontSize',12);
set(gca,'FontSize',12);
title('Measurements \tau_c(s) in seconds ','FontSize',12);

subplot(1,3,2)
hh=scatter(slon(IJ),slat(IJ),100,2*sqrt(G*m2),'filled'); hold on;
plot(slon(IJ),slat(IJ),'ko','MarkerSize',10);
colormap(jet); colorbar; caxis([0.15 0.50])
plot(lone,late,'w^','MarkerSize',15,'MarkerFaceColor','r');
axis([-117.0 -115.9 33.0 34.1])
xlabel('Longitude','FontSize',12);
ylabel('Latitude','FontSize',12);
title('2nd Moments Fit, \tau_c(s) in seconds','FontSize',12);
set(gca,'FontSize',12);

subplot(1,3,3)
[xq, yq] = meshgrid(min(slon(IJ)):.001:max(slon(IJ)), min(slat(IJ)):.001:max(slat(IJ)));
vq = griddata(slon(IJ), slat(IJ), 2*sqrt(G*m2), xq, yq);
contour(xq,yq,vq);
colormap(jet); colorbar; caxis([0.15 0.50])
xlabel('Longitude')
ylabel('Latitude')

%plot(slon(IJ), slat(IJ), 2*sqrt(G*m2), 'filled');hold on;
%colormap(jet);colorbar

%[c,h] = contour(xq,yq,vq);
%clabel(c,h)

if mode_run.debug_plot
    k = waitforbuttonpress;
    close
else
    set(gcf, 'visible', 'off')
    set(gcf, 'pos', [10 10 900 500])
end

if ~mode_run.no_figure
    saveas(gcf, sprintf('%s/MS%d_EGF%d_FITresult.png', image_dir, msorid, orid))
end

