
variables = {'tos','tas','pr','psl','monmaxpr','monmaxtasmax','monmintasmin','zmta'};
season = 1:12;

member = {'1A','1B','1C','1D','1E','1F','1G','1H','1J'};

for k = 1:length(variables)
    variable = variables{k};

    submission_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/submissions-Tier1-standardized-estimates');
    emean_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/ensmeans-Tier1');
    ref_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/Evaluation-Tier1');

    %member = {'1B','1D','1E','1G','1J'};
    %member = {'1B','1E','1G','1J'};
    startyear = 1980;
    endyear = 2022;

    submission_short_names = {'RegGMST','4th-Order-Polynomial','10yr-Lowpass','LFCA','LFCA-2', ...
'MF-LFCA','MF-LFCA-2','LIMnMCA','ICA-lowpass','LIMopt','LIMopt-filter', ...
'Colored-LIMnMCA','DMDc','GPCA','GPCA-DA','RegGMST-LENSem','MLR-Forcing', ...
'SNMP-OF','AllFinger','MonthFinger','3DUNet-Fingerprinters','EOF-SLR','LDM-SLR',...
'Anchor-OPLS','UNet3D-LOCEAN','TrainingEM','RandomForest','EncoderDecoder','EnsFMP','ANN-Fingerprinters'};

    switch variable
        case 'monmaxtasmax'
            varnam = 'tasmax';
        case 'monmaxpr'
            varnam = 'pr';
        case 'monmintasmin'
            varnam = 'tasmin';
        case 'zmta'
            varnam = 'ta';
        otherwise
            varnam = variable;
    end

    clear rmse_ref rmse_all corr_ref corr_all std_ref std_all std_emean

    n1 = startyear - 1950 + 1;
    n2 = endyear - 1950 + 1;

    clear fields_seasonal

    for j = 1:length(member)
        if strcmp(variable,'zmta')
            %simplicity_order = [nan 2 1 12 13 17 18 14 nan 10 15 16 4 8 9 nan 3 19 nan nan 20 5 11 nan 22 21 nan 6 7 nan];
            % with RegGMST and RegGMST-LENSem switched (until data files are updated)
            simplicity_order = [nan 2 1 12 13 17 18 14 10 15 16 4 nan 8 9 nan 19 20 3 nan nan 5 11 nan 22 21 nan 6 7 nan];
        else
            %simplicity_order = [26 3 1 16 17 21 22 18 8 14 19 20 7 12 13 27 5 23 24 2 28 9 15 6 30 29 25 10 11 4];
            % with RegGMST and RegGMST-LENSem manually switched (until data files are updated)
            simplicity_order = [27 3 1 16 17 21 22 18 14 19 20 7 8 12 13 26 23 28 5 24 2 9 15 6 30 29 25 10 11 4];
        end

        submission_file = select_files(submission_files,[variable,'_']);
        submission_file = select_files(submission_file,member{j});
        emean_file = select_files(emean_files,[variable,'.']);
        if strcmp(member{j},'1I')
            disp('Member 1I is observations, no ensemble mean found.')
            emean_file = select_files(emean_file,'1A');
        else
            emean_file = select_files(emean_file,member{j});
        end
        ref_file = select_files(ref_files,[variable,'_']);
        ref_file = select_files(ref_file,member{j});
        if strcmp(variable,'pr')
            submission_file = select_files(submission_file,member{j},'monmax');
            emean_file = select_files(emean_file,member{j},'monmax');
            ref_file = select_files(ref_file,member{j},'monmax');
        end
        if strcmp(variable,'zmta')
            [fields,lat,plev,time] = get_avg_field_nd(submission_file,{'forced_component','lat','plev','time'});
        else
            [fields,lat,lon,time] = get_avg_field_nd(submission_file,{'forced_component','lat','lon','time'});
        end
        try
            field_emean = get_avg_field_nd(emean_file,'arr_EM');
        catch
            field_emean = get_avg_field_nd(emean_file,'tos'); % not necessary if I get this file (1E tos) from Adam
        end
        field_emean = squeeze(field_emean);
        field_ref = get_avg_field_nd(ref_file,varnam);

        fields(abs(fields)>1e10) = nan;
        fields(fields==0) = nan;
        field_emean(abs(field_emean)>1e10) = nan;
        field_emean(field_emean==0) = nan;
        field_ref(abs(field_ref)>1e10) = nan;
        field_ref(field_ref==0) = nan;
        if strcmp(variable,'zmta')
            tmp = fields;
            for i = 1:length(simplicity_order)
                if isnan(simplicity_order(i))
                    tmp(:,:,:,i) = nan;
                else
                    tmp(:,:,:,i) = fields(:,:,:,simplicity_order(i));
                end
            end
            fields = tmp; clear tmp
        else
            fields = fields(:,:,:,simplicity_order);
        end
        s = size(fields);

        nyr = s(3)/12;
        months = repmat(1:12,[1 nyr]);
        years = floor(1950+1/24:1/12:2022.99);

        clear clim_ref

        for n = 1:12
            clim_ref(:,:,n) = mean(field_ref(:,:,n:12:end),3);
            field_ref(:,:,n:12:end) = field_ref(:,:,n:12:end) - mean(field_ref(:,:,n:12:end),3);
            field_emean(:,:,n:12:end) = field_emean(:,:,n:12:end) - mean(field_emean(:,:,n:12:end),3);
            fields(:,:,n:12:end,:) = fields(:,:,n:12:end,:) - mean(fields(:,:,n:12:end,:),3);
        end

        if contains(variable,'max')
            field_ref = monthly_to_seasonal(months,years,field_ref+repmat(clim_ref,[1 1 s(3)/12]),season,3,'max');
            field_emean = monthly_to_seasonal(months,years,field_emean+repmat(clim_ref,[1 1 s(3)/12]),season,3,'max');
        elseif contains(variable,'min')
            field_ref = monthly_to_seasonal(months,years,field_ref+repmat(clim_ref,[1 1 s(3)/12]),season,3,'min');
            field_emean = monthly_to_seasonal(months,years,field_emean+repmat(clim_ref,[1 1 s(3)/12]),season,3,'min');
        else
            field_ref = monthly_to_seasonal(months,years,field_ref,season,3);
            field_emean = monthly_to_seasonal(months,years,field_emean,season,3);
        end
        for i = 1:s(4)
            if contains(variable,'max')
                fields_seasonal(:,:,:,i) = monthly_to_seasonal(months,years,fields(:,:,:,i)+repmat(clim_ref,[1 1 s(3)/12]),season,3,'max');
            elseif contains(variable,'min')
                fields_seasonal(:,:,:,i) = monthly_to_seasonal(months,years,fields(:,:,:,i)+repmat(clim_ref,[1 1 s(3)/12]),season,3,'min');
            else
                fields_seasonal(:,:,:,i) = monthly_to_seasonal(months,years,fields(:,:,:,i),season,3);
            end
        end

        [~,trend_ref] = detrend(field_ref(:,:,n1:n2),3,1);
        [~,trend_emean] = detrend(field_emean(:,:,n1:n2),3,1);
        [~,trends] = detrend(fields_seasonal(:,:,n1:n2,:),3,1);
        trends = squeeze(trends);

        if strcmp(variable,'tos') || strcmp(variable,'zmta')
            % make sure NaNs in all the same land places
            mask = trends(:,:,2)./trends(:,:,2);
            mask_ref = trend_ref./trend_ref;
            mask_emean = trend_emean./trend_emean;
            trend_ref = trend_ref.*mask.*mask_ref.*mask_emean;
            trends = trends.*mask.*mask_ref.*mask_emean;
            trend_emean = trend_emean.*mask.*mask_ref.*mask_emean;
        end

        for i = 1:s(4)
            if strcmp(variable,'zmta')
                rmse_ref(j) = sqrt(global_mean(plev,lat,((trend_emean-trend_ref).^2)'));
                rmse_all(i,j) = sqrt(global_mean(plev,lat,((trend_emean-trends(:,:,i)).^2)'));
                std_emean(j) = sqrt(global_mean(plev,lat,(trend_emean.^2)'));
                std_ref(j) = sqrt(global_mean(plev,lat,(trend_ref.^2)'));
                std_all(i,j) = sqrt(global_mean(plev,lat,(trends(:,:,i).^2)'));
                corr_ref(j) = global_mean(plev,lat,(trend_emean.*trend_ref)')./sqrt(global_mean(plev,lat,(trend_emean.^2)').*global_mean(plev,lat,(trend_ref.^2)'));
                corr_all(i,j) = global_mean(plev,lat,(trend_emean.*trends(:,:,i))')./sqrt(global_mean(plev,lat,(trend_emean.^2)').*global_mean(plev,lat,(trends(:,:,i).^2)'));
            else
                rmse_ref(j) = sqrt(global_mean(lon,lat,(trend_emean-trend_ref).^2));
                rmse_all(i,j) = sqrt(global_mean(lon,lat,(trend_emean-trends(:,:,i)).^2));
                std_emean(j) = sqrt(global_mean(lon,lat,trend_emean.^2));
                std_ref(j) = sqrt(global_mean(lon,lat,trend_ref.^2));
                std_all(i,j) = sqrt(global_mean(lon,lat,trends(:,:,i).^2));
                corr_ref(j) = global_mean(lon,lat,trend_emean.*trend_ref)./sqrt(global_mean(lon,lat,trend_emean.^2).*global_mean(lon,lat,trend_ref.^2));
                corr_all(i,j) = global_mean(lon,lat,trend_emean.*trends(:,:,i))./sqrt(global_mean(lon,lat,trend_emean.^2).*global_mean(lon,lat,trends(:,:,i).^2));
            end
        end

    end

    std_emean_all(:,k) = std_emean;
    rmse_ref_all(:,k) = rmse_ref;
    corr_ref_all(:,k) = corr_ref;
    rmse_all_all(:,:,k) = rmse_all;
    corr_all_all(:,:,k) = corr_all;

    % normalized before averaging
    std_ref_combined(k) = sqrt(mean((std_ref./std_emean).^2));
    std_all_combined(:,k) = sqrt(mean((std_all./std_emean).^2,2));

    rmse_ref_combined(k) = sqrt(mean((rmse_ref./std_emean).^2));
    rmse_all_combined(:,k) = sqrt(mean((rmse_all./std_emean).^2,2));

    std_emean_combined(k) = sqrt(mean(std_emean.^2));

    % not used in Taylor diagram
    corr_ref_combined(k) = sqrt(mean(corr_ref.^2));
    corr_all_combined(:,k) = sqrt(mean(corr_all.^2,2));

    k
end

%%

variable_labels = {'SST','T2m','PR','SLP','Rx1day','TXx','TNn','zmTa'};

nme = size(rmse_all_all,1);

unseen = [2 4 5 7 9]; 
nmo = length(unseen);
nv = length(variables);
rmse_all_c = squeeze(sqrt(nanmean((rmse_all_all(:,unseen,:)./repmat(reshape(std_emean_all(unseen,:),[1 nmo nv]),[nme 1 1])).^2,2)));
rmse_ref_c = sqrt(nanmean((rmse_ref_all(unseen,:)./std_emean_all(unseen,:)).^2));
% corr_all_c = squeeze(sqrt(nanmean(corr_all_all(:,unseen,:).^2,2)));
% corr_ref_c = sqrt(nanmean(corr_ref_all(unseen,:).^2));
corr_all_c = squeeze(nanmean(corr_all_all(:,unseen,:),2));
corr_ref_c = nanmean(corr_ref_all(unseen,:));

%to_plot = 1-(rmse_all_combined./rmse_ref_combined).^2;
%to_plot = 1-(rmse_all_c./rmse_ref_c).^2;
to_plot = [(1-rmse_ref_c); 1-rmse_all_c];
stipple = to_plot>(1-rmse_ref_c);
to_plot(:,end+1) = nan;
to_plot(end+1,:) = nan;
plot_field_pixels(1:length(submission_short_names)+2,1:9,to_plot,linspace(-1,1,25)) % -1,1,21
plot_stipple(1.5:length(submission_short_names)+1.5,1.5:8.5,stipple,6)
set(gca,'color',[0.6 0.6 0.6])
pretty_figure(700,250,'none','none','none','none',16)
title('1 - RMSE (normalized)','fontsize',16)
set(gca,'ytick',1.5:1:8.5); set(gca,'yticklabel',variable_labels);
set(gca,'ydir','reverse'); set(gca,'box','off')
set(gca,'xtick',[1.5:1:11.5, 13.5:2:length(submission_short_names)+1.5]); set(gca,'xticklabel',[0:1:10, 12:2:length(submission_short_names)]);
xtickangle(0)

%to_plot = corr_all_combined-corr_ref_combined;
%to_plot = corr_all_c-corr_ref_c;
to_plot = [corr_ref_c; corr_all_c];
stipple = to_plot>corr_ref_c;
to_plot(:,end+1) = nan;
to_plot(end+1,:) = nan;
plot_field_pixels(1:length(submission_short_names)+2,1:9,to_plot.^2,linspace(0,1,26)) % -0.2,0.2,21
hc = colorbar; set(hc,'ytick',[0.2 0.4 0.6 0.7 0.8 0.9 1].^2); set(hc,'yticklabel',[0.2 0.4 0.6 0.7 0.8 0.9 1])
plot_stipple(1.5:length(submission_short_names)+1.5,1.5:8.5,stipple,6)
set(gca,'color',[0.6 0.6 0.6])
pretty_figure(700,250,'none','none','none','none',16)
title('Pattern Correlation','fontsize',16)
set(gca,'ytick',1.5:1:8.5); set(gca,'yticklabel',variable_labels);
set(gca,'ydir','reverse'); set(gca,'box','off')
set(gca,'xtick',[1.5:1:11.5, 13.5:2:length(submission_short_names)+1.5]); set(gca,'xticklabel',[0:1:10, 12:2:length(submission_short_names)]);
xtickangle(0)

%%

unseen = [2 4 5 7 9];
seen = [1 3 6 8];
is = [3 6 8 1 2 4 5 7 9];
models = {'CanESM5','CESM2','MIROC6','MPI-ESM1-2-LR','ACCESS-ESM1-5','EC-Earth3','GFDL-SPEAR-MED','IPSL-CM6A-LR','NorCPM1'};

for i = 1:length(variables)
    rmse = rmse_all_all(:,is,i)./std_emean_all(is,i)';
    rmse_ref = rmse_ref_all(is,i)'./std_emean_all(is,i)';
    %to_plot = 1-(rmse./rmse_ref).^2;
    to_plot = [1-rmse_ref; 1-rmse];
    stipple = to_plot>(1-rmse_ref);
    to_plot(:,end+1) = nan;
    to_plot(end+1,:) = nan;
    plot_field_pixels(1:length(submission_short_names)+2,1:10,to_plot,linspace(-1,1,25))
    plot_stipple(1.5:length(submission_short_names)+1.5,1.5:9.5,stipple,6)
    set(gca,'color',[0.6 0.6 0.6])
    pretty_figure(700,250,'none','none','none','none',16)
    title([variable_labels{i},', 1 - RMSE (normalized)'],'fontsize',16)
    set(gca,'ytick',1.5:1:9.5); set(gca,'yticklabel',models);
    set(gca,'ydir','reverse'); set(gca,'box','off')
    set(gca,'xtick',[1.5:1:11.5, 13.5:2:length(submission_short_names)+1.5]); set(gca,'xticklabel',[0:1:10, 12:2:length(submission_short_names)]);
    xtickangle(0)
    hold on; plot([1 32],[5 5],'k')

    %to_plot = corr_all_all(:,is,i)-corr_ref_all(is,i)';
    to_plot = [corr_ref_all(is,i)'; corr_all_all(:,is,i)];
    stipple = to_plot>corr_ref_all(is,i)';
    to_plot(:,end+1) = nan;
    to_plot(end+1,:) = nan;
    if sum(i == [3 4 5])==1
        plot_field_pixels(1:length(submission_short_names)+2,1:10,to_plot,linspace(-0.05,1,22))
    elseif i == 8
        plot_field_pixels(1:length(submission_short_names)+2,1:10,to_plot,linspace(0.79,1,22))
    else
        plot_field_pixels(1:length(submission_short_names)+2,1:10,to_plot,linspace(0.575,1,22))
    end
    plot_stipple(1.5:length(submission_short_names)+1.5,1.5:9.5,stipple,6)
    set(gca,'color',[0.6 0.6 0.6])
    pretty_figure(700,250,'none','none','none','none',16)
    title([variable_labels{i},', Pattern Correlation'],'fontsize',16)
    set(gca,'ytick',1.5:1:9.5); set(gca,'yticklabel',models);
    set(gca,'ydir','reverse'); set(gca,'box','off')
    set(gca,'xtick',[1.5:1:11.5, 13.5:2:length(submission_short_names)+1.5]); set(gca,'xticklabel',[0:1:10, 12:2:length(submission_short_names)]);
    xtickangle(0)
    hold on; plot([1 32],[5 5],'k')
end

