
% Used to make Figs. A1-A3 of Wills et al. 2025
% to be run after ForceSMIP_trend_skill.m

index = 1:30;
estimates = trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2))));
estimates_anom = estimates - mean(estimates,3);
estimates_anom(isnan(estimates_anom)) = 0;

[EOFs,PCs,pvar] = EOF_analysis(lon,lat,estimates_anom);

pvar(1)
pvar(2)

switch variable
    case 'tos'
        EOFs = -EOFs; PCs = -PCs;
        plot_field_robinson(lon,lat,std(trends,0,3),linspace(-0.8,0.8,25));
        plot_field_robinson(lon,lat,EOFs(:,:,1),linspace(-0.8,0.8,25));
        plot_field_robinson(lon,lat,EOFs(:,:,2),linspace(-0.8,0.8,25));
    case 'psl'
        EOFs = -EOFs; PCs = -PCs;
        plot_field_robinson(lon,lat,std(trends,0,3),linspace(-200,200,25));
        plot_field_robinson(lon,lat,EOFs(:,:,1),linspace(-200,200,25));
        plot_field_robinson(lon,lat,EOFs(:,:,2),linspace(-200,200,25));
    case 'pr'
        plot_field_robinson(lon,lat,std(trends,0,3),linspace(-1.6,1.6,25),[0 180 0],[-90 90],[0 360],'none','none','precip');
        plot_field_robinson(lon,lat,EOFs(:,:,1),linspace(-1.6,1.6,25),[0 180 0],[-90 90],[0 360],'none','none','precip');
        plot_field_robinson(lon,lat,EOFs(:,:,2),linspace(-1.6,1.6,25),[0 180 0],[-90 90],[0 360],'none','none','precip');
end

labels = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
labels = labels(index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2))));

i = 1;
while str2num(labels{i}) < 18
    i = i + 1;
end
icomplex = i;

dx = 0.05;
dy = 0.05;

figure; plot(PCs(1:(icomplex-1),1),PCs(1:(icomplex-1),2),'ko','markerfacecolor','k');
hold on; plot(PCs(icomplex:end,1),PCs(icomplex:end,2),'kd','markerfacecolor','k');
text(PCs(:,1)+dx,PCs(:,2)+dy,labels,'fontsize',12)
%pretty_figure(500,400,'PC1','PC2','none','none',16);
set(gca,'xgrid','on'); set(gca,'ygrid','on')
