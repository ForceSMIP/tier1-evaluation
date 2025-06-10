
% Used to make Figs. 11 and S3 of Wills et al. 2025

variables = {'tos','tas','pr','psl','monmaxpr','monmaxtasmax','monmintasmin','zmta'};

index = 1:30;

for k = 1:length(variables)
    variable = variables{k};
    ForceSMIP_trend_skill % make sure first line (specifying variable) is commented out
    switch variable
        case 'tos'
            ctrs = linspace(-2,2,25);
        case 'tas'
            ctrs = linspace(-3,3,25);
        case 'psl'
            ctrs = linspace(-240,240,25);
        case 'pr'
            ctrs = linspace(-1.6,1.6,25);
        case 'monmaxtasmax'
            ctrs = linspace(-3,3,25);
        case 'monmintasmin'
            ctrs = linspace(-3,3,25);
        case 'monmaxpr'
            ctrs = linspace(-20,20,25);
        case 'zmta'
            ctrs = linspace(-2,2,25);
    end
      
    if strcmp(variable,'zmta')
%         plot_field_div(lat,plev./100,trends(:,:,26),ctrs);
%         ylabel('Pressure (hPa)')
%         set(gca,'xtick',-60:30:60); set(gca,'xticklabel',{'60°S','30°S','EQ','30°N','60°N'})
%         set(gca,'color',[0.8 0.8 0.8])
%         set(gca,'ydir','reverse')
        plot_field_div(lat,plev./100,mean(trends(:,:,index(RMSs(3:end)<RMSs(2)&STDs(3:end)>0.5&CORs(3:end)>0.3)),3),ctrs);
        ylabel('Pressure (hPa)')
        set(gca,'xtick',-60:30:60); set(gca,'xticklabel',{'60°S','30°S','EQ','30°N','60°N'})
        set(gca,'color',[0.8 0.8 0.8])
        set(gca,'ydir','reverse')
        plot_field_div(lat,plev./100,trend_ref-mean(trends(:,:,index(RMSs(3:end)<RMSs(2)&STDs(3:end)>0.5&CORs(3:end)>0.3)),3),ctrs);
        ylabel('Pressure (hPa)')
        set(gca,'xtick',-60:30:60); set(gca,'xticklabel',{'60°S','30°S','EQ','30°N','60°N'})
        set(gca,'color',[0.8 0.8 0.8])
        set(gca,'ydir','reverse')
    elseif strfind(variable,'pr')
        %plot_field_robinson(lon,lat,trends(:,:,26),ctrs,[0 180 0],[-90 90],[0 360],'none','none','precip');
        plot_field_robinson(lon,lat,mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs,[0 180 0],[-90 90],[0 360],'none','none','precip');
        plot_field_robinson(lon,lat,trend_ref-mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs,[0 180 0],[-90 90],[0 360],'none','none','precip');
    else
        %plot_field_robinson(lon,lat,trends(:,:,26),ctrs);
        plot_field_robinson(lon,lat,mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs);
        plot_field_robinson(lon,lat,trend_ref-mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs);
    end

    k
end