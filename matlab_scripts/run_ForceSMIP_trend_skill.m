
% Used to make Figs. 11 and S3 of Wills et al. 2025

variables = {'tos','tas','pr','psl','monmaxpr','monmaxtasmax','monmintasmin','zmta'};

index = 1:30;

for k = 1:length(variables)
    variable = variables{k};
    ForceSMIP_trend_skill % make sure first line (specifying variable) is commented out
    switch variable
        case 'tos'
            ctrs = linspace(-2,2,25);
            tos_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            tos_trainingEM = fields_seasonal(:,:,:,26);
            tos_raw = field_ref;
        case 'tas'
            ctrs = linspace(-3,3,25);
            tas_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            tas_trainingEM = fields_seasonal(:,:,:,26);
            tas_raw = field_ref;
        case 'psl'
            ctrs = linspace(-240,240,25);
            psl_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            psl_trainingEM = fields_seasonal(:,:,:,26);
            psl_raw = field_ref;
        case 'pr'
            ctrs = linspace(-1.6,1.6,25);
            pr_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            pr_trainingEM = fields_seasonal(:,:,:,26);
            pr_raw = field_ref;
        case 'monmaxtasmax'
            ctrs = linspace(-3,3,25);
            monmaxtasmax_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            monmaxtasmax_trainingEM = fields_seasonal(:,:,:,26);
            monmaxtasmax_raw = field_ref;
        case 'monmintasmin'
            ctrs = linspace(-3,3,25);
            monmintasmin_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            monmintasmin_trainingEM = fields_seasonal(:,:,:,26);
            monmintasmin_raw = field_ref;
        case 'monmaxpr'
            ctrs = linspace(-20,20,25);
            monmaxpr_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            monmaxpr_trainingEM = fields_seasonal(:,:,:,26);
            monmaxpr_raw = field_ref;
        case 'zmta'
            ctrs = linspace(-2,2,25);
            zmta_method_mean = mean(fields_seasonal(:,:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),4);
            zmta_trainingEM = fields_seasonal(:,:,:,26);
            zmta_raw = field_ref;
    end
      
    if strcmp(variable,'zmta')
%         plot_field_div(lat,plev./100,trends(:,:,26),ctrs);
%         ylabel('Pressure (hPa)')
%         set(gca,'xtick',-60:30:60); set(gca,'xticklabel',{'60°S','30°S','EQ','30°N','60°N'})
%         set(gca,'color',[0.8 0.8 0.8])
%         set(gca,'ydir','reverse')
        plot_field_div(lat,plev./100,mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs);
        ylabel('Pressure (hPa)')
        set(gca,'xtick',-60:30:60); set(gca,'xticklabel',{'60°S','30°S','EQ','30°N','60°N'})
        set(gca,'color',[0.8 0.8 0.8])
        set(gca,'ydir','reverse')
        plot_field_div(lat,plev./100,trend_ref-mean(trends(:,:,index((1-CORs(3:end)/CORs(2))<(1-RMSs(3:end)/RMSs(2)))),3),ctrs);
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

year = 1950:2022;
ncbulkexport_multidimensional('ForceSMIP_Tier1_OBS_estimate.nc',lon,'lon',lat,'lat',year,'year',{tos_raw,tos_trainingEM,tos_method_mean,tas_raw,tas_trainingEM,tas_method_mean, ...
    psl_raw,psl_trainingEM,psl_method_mean,pr_raw,pr_trainingEM,pr_method_mean,monmaxtasmax_raw,monmaxtasmax_trainingEM,monmaxtasmax_method_mean, ...
    monmintasmin_raw,monmintasmin_trainingEM,monmintasmin_method_mean,monmaxpr_raw,monmaxpr_trainingEM,monmaxpr_method_mean}, ...
    {'tos_raw','tos_trainingEM','tos_method_mean','tas_raw','tas_trainingEM','tas_method_mean', ...
    'psl_raw','psl_trainingEM','psl_method_mean','pr_raw','pr_trainingEM','pr_method_mean','monmaxtasmax_raw','monmaxtasmax_trainingEM','monmaxtasmax_method_mean', ...
    'monmintasmin_raw','monmintasmin_trainingEM','monmintasmin_method_mean','monmaxpr_raw','monmaxpr_trainingEM','monmaxpr_method_mean'});

