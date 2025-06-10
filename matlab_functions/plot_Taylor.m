function [h] = plot_Taylor(STDs,STD_ref,CORs,labels,xmin,xmax,ymax,cmax,cticks)

STDs = STDs./STD_ref;

alpha = acos(CORs);
x = STDs.*cos(alpha);
y = STDs.*sin(alpha);
labels = labels(~isnan(x));
x = x(~isnan(x));
y = y(~isnan(y));

if nargin < 5
    ymax = ceil(max(y*10))/10;
    xmax = ceil(max(x*10))/10;
    xmin = min(1,floor(min(x*10))/10);
end

xi = xmin:0.01:xmax;
yi = 0:0.01:ymax;

xi_hi = xmin:0.001:xmax;
yi_hi = 0:0.001:ymax;

[Y,X] = meshgrid(yi,xi);
ci = sqrt((X-1).^2+Y.^2);
ci2 = sqrt(X.^2+Y.^2);

if nargin < 8
    cmax = ceil(max(max(ci))*10)/10;
    cticks = 21;
end

plot_field(xi,yi,ci,linspace(0,cmax,cticks));
set(gca,'ylim',[0 ymax])
set(gca,'xlim',[xmin xmax])

yticks = []; yticklabels = [];
for i = [-0.99 -0.95 -0.9:0.1:0.9, 0.95 0.99]
    ycorr = xi.*tan(acos(i));
    hold on; plot(xi,ycorr,'k')
    if xmax*tan(acos(i))<ymax & xmax*tan(acos(i))>0
        yticks = [yticks xmax*tan(acos(i))];
        yticklabels = [yticklabels i];
    end
end
hold on; plot([0 0],[0 ymax],'k','linewidth',1.5)

if xmax<3
    hold on; contour(xi,yi,ci2',0.2:0.2:10,'k')
else
    hold on; contour(xi,yi,ci2',0.5:0.5:10,'k')
end
hold on; contour(xi,yi,ci2',[1 1],'k','linewidth',1.5)

%hold on; plot(xi,xi.*tan(alpha(1)),':','color','w','linewidth',1.5);
%hold on; plot(xi,sqrt((x(1)-1).^2+y(1).^2-(xi-1).^2),':','color','w','linewidth',1.5);

fred_rmse = 1-sqrt((xi_hi'-1).^2+yi_hi.^2)./sqrt((x(1)-1).^2+y(1).^2);
fred_pcorr = 1-(xi_hi'./sqrt(xi_hi'.^2+yi_hi.^2))./(x(1)./sqrt(x(1)^2+y(1)^2));
hold on; contour(xi_hi',yi_hi,(fred_rmse>fred_pcorr)',[0.5 0.5],':','color','w','linewidth',1.5);

% here the information about the type of tier 1 methods is hard-coded
hold on; plot(x(2:18),y(2:18),'ko','markerfacecolor','k','markersize',6);
hold on; plot(x(19:end),y(19:end),'kd','markerfacecolor','k','markersize',6);

if nargin > 3
    if xmax<2
        dx = xmax/160;
        dy = ymax/160;
    else
        dx = xmax/80;
        dy = ymax/80;
    end
    text(x+dx,y+dy,labels,'fontsize',13)
end

hold on; plot(x(1),y(1),'kpentagram','markerfacecolor','w','markersize',11,'linewidth',1); %[0.9 0.9 0.9]

%pretty_figure(625,550,'none','none','none','none',16);
pretty_figure(620,560,'none','none',-0.2:0.2:2.4,fliplr(yticks),16,-0.2:0.2:2.4,fliplr(yticklabels));
set(gca,'YAxisLocation','right')
set(gca,'box','off')
hold on; line([xmin xmin],[0 ymax],'Color','k')
hold on; line([xmin xmax],[ymax ymax],'Color','k')
colorbar off; colorbar('southoutside')