
% Used to make Figs. 1j, 2j, 3, 4, and S2 of Wills et al. 2025
% The other panels of Figs. 1 and 2 can be made by pausing at l. 161 and
% plotting trend_ref, trends(:,:,i), and trends(:,:,i)-trend_ref
% Figs. 8-10 can be made by plotting trend_ref, trends(:,:,i), and trends(:,:,i)-trend_ref

variable = 'zmta'; % 'tos','tas','pr','psl','monmaxpr','monmaxtasmax','monmintasmin','zmta'
season = 1:12; % for seaonal mean (1:12 makes annual mean)

clear fields_seasonal

submission_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/submissions-Tier1-standardized-estimates');
emean_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/ensmeans-Tier1');
ref_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/Evaluation-Tier1');

member = {'1B','1D','1E','1G','1J'}; % unseen models '1B','1D','1E','1G','1J' in ForceSMIP Tier 1
obs = {'1I'};
startyear = 1980;
endyear = 2022;

submission_short_names = {'RegGMST','4th-Order-Polynomial','10yr-Lowpass','LFCA','LFCA-2', ...
'MF-LFCA','MF-LFCA-2','LIMnMCA','ICA-lowpass','LIMopt','LIMopt-filter', ...
'Colored-LIMnMCA','DMDc','GPCA','GPCA-DA','RegGMST-LENSem','MLR-Forcing', ...
'SNMP-OF','AllFinger','MonthFinger','3DUNet-Fingerprinters','EOF-SLR','LDM-SLR',...
'Anchor-OPLS','UNet3D-LOCEAN','TrainingEM','RandomForest','EncoderDecoder','EnsFMP','ANN-Fingerprinters'};

switch variable % some variables have different variable names than file names
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

n1 = startyear - 1950 + 1;
n2 = endyear - 1950 + 1;

clear fields_seasonal

rmse_ref = zeros(1,length(member)); std_emean = rmse_ref; std_ref = rmse_ref; corr_ref = rmse_ref;
rmse_all = zeros(s(4),length(member)); std_all = rmse_all; corr_all = rmse_all;

for j = 1:length(member)
    
    % specifying simplicity ordering used in Wills et al. 2025
    if strcmp(variable,'zmta')
        %simplicity_order = [nan 2 1 12 13 17 18 14 nan 10 15 16 4 8 9 nan 3 19 nan nan 20 5 11 nan 22 21 nan 6 7 nan];
        % with RegGMST and RegGMST-LENSem switched (until data files are updated)
        simplicity_order = [nan 2 1 12 13 17 18 14 10 15 16 4 nan 8 9 nan 19 20 3 nan nan 5 11 nan 22 21 nan 6 7 nan];
    else
        %simplicity_order = [26 3 1 16 17 21 22 18 8 14 19 20 7 12 13 27 5 23 24 2 28 9 15 6 30 29 25 10 11 4];
        % with RegGMST and RegGMST-LENSem manually switched (until data files are updated)
        simplicity_order = [27 3 1 16 17 21 22 18 14 19 20 7 8 12 13 26 23 28 5 24 2 9 15 6 30 29 25 10 11 4];
    end

    % find all filenames
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
    if strcmp(variable,'pr') % remove monmaxpr files from file lists
        submission_file = select_files(submission_file,member{j},'monmax');
        emean_file = select_files(emean_file,member{j},'monmax');
        ref_file = select_files(ref_file,member{j},'monmax');
    end

    % load data, change missing data into NaN, apply simplicity ordering
    if strcmp(variable,'zmta')
        [fields,lat,plev,time] = get_avg_field_nd(submission_file,{'forced_component','lat','plev','time'});
    else
        [fields,lat,lon,time] = get_avg_field_nd(submission_file,{'forced_component','lat','lon','time'});
    end
    field_emean = get_avg_field_nd(emean_file,'arr_EM');
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
    years = floor(1950:1/12:2022.99);

    clear clim_ref

    % compute seasonal climatologies and anomalies
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

    % compute linear trends (in units per period length) 
    [~,trend_ref] = detrend(field_ref(:,:,n1:n2),3,1);
    [~,trend_emean] = detrend(field_emean(:,:,n1:n2),3,1);
    [~,trends] = detrend(fields_seasonal(:,:,n1:n2,:),3,1);
    trends = squeeze(trends);

    % apply uniform NaNs to all fields (in case of any differences)
    if strcmp(variable,'tos') || strcmp(variable,'zmta')
        % make sure NaNs in all the same land places
        mask = trends(:,:,2)./trends(:,:,2);
        mask_ref = trend_ref./trend_ref;
        mask_emean = trend_emean./trend_emean;
        trend_ref = trend_ref.*mask.*mask_ref.*mask_emean;
        trends = trends.*mask.*mask_ref.*mask_emean;
        trend_emean = trend_emean.*mask.*mask_ref.*mask_emean;
    end

    % compute skill metrics
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

%% Process observations

% find files
submission_file = select_files(submission_files,[variable,'_']);
submission_file = select_files(submission_file,obs);
ref_file = select_files(ref_files,[variable,'_']);
ref_file = select_files(ref_file,obs);
if strcmp(variable,'pr')
    submission_file = select_files(submission_file,obs,'monmax');
    ref_file = select_files(ref_file,obs,'monmax');
end
% load data, replace missing values with NaN, apply simplicity ordering
if strcmp(variable,'zmta')
    [fields,lat,plev,time] = get_avg_field_nd(submission_file,{'forced_component','lat','plev','time'});
else
    [fields,lat,lon,time] = get_avg_field_nd(submission_file,{'forced_component','lat','lon','time'});
end
field_ref = get_avg_field_nd(ref_file,varnam);
fields(abs(fields)>1e10) = nan;
fields(fields==0) = nan;
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

% compute seasonal climatologies and anomalies
nyr = s(3)/12;
months = repmat(1:12,[1 nyr]);
years = floor(1950+1/24:1/12:2022.99);

for n = 1:12
    clim_ref(:,:,n) = mean(field_ref(:,:,n:12:end),3);
    field_ref(:,:,n:12:end) = field_ref(:,:,n:12:end) - mean(field_ref(:,:,n:12:end),3);
    fields(:,:,n:12:end,:) = fields(:,:,n:12:end,:) - mean(fields(:,:,n:12:end,:),3);
end

if contains(variable,'max')
    field_ref = monthly_to_seasonal(months,years,field_ref+repmat(clim_ref,[1 1 s(3)/12]),season,3,'max');
elseif contains(variable,'min')
    field_ref = monthly_to_seasonal(months,years,field_ref+repmat(clim_ref,[1 1 s(3)/12]),season,3,'min');
else
    field_ref = monthly_to_seasonal(months,years,field_ref,season,3);
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

% compute linear trends (in units per period length) 
[~,trend_ref] = detrend(field_ref(:,:,n1+1:n2),3,1);
[~,trends] = detrend(fields_seasonal(:,:,n1+1:n2,:),3,1);
trends = squeeze(trends);

% apply uniform NaNs to all fields (in case of any differences)
if strcmp(variable,'tos') || strcmp(variable,'zmta')
    % make sure NaNs in all the same land places
    mask = trends(:,:,2)./trends(:,:,2);
    mask_ref = trend_ref./trend_ref;
    mask_emean = trend_emean./trend_emean;
    trend_ref = trend_ref.*mask.*mask_ref.*mask_emean;
    trends = trends.*mask.*mask_ref.*mask_emean;
    trend_emean = trend_emean.*mask.*mask_ref.*mask_emean;
end

%% averaging skill metrics

% normalized before averaging
std_ref = sqrt(mean((std_ref./std_emean).^2));
std_all = sqrt(mean((std_all./std_emean).^2,2));

rmse_ref = sqrt(mean((rmse_ref./std_emean).^2));
rmse_all = sqrt(mean((rmse_all./std_emean).^2,2));

std_emean = sqrt(mean(std_emean.^2));

% not used in Taylor diagram
corr_ref = sqrt(mean(corr_ref.^2));
corr_all = sqrt(mean(corr_all.^2,2));

%% Plot Taylor diagram

STDs = [1, std_ref, std_all'];
%CORs = [1, corr_ref, corr_all'];
%RMSs(1) = 0;
%RMSs(2:s(4)+2) = sqrt(STDs(2:s(4)+2).^2 + STDs(1)^2 - 2*STDs(2:s(4)+2)*STDs(1).*CORs(2:s(4)+2));
RMSs = [0 rmse_ref, rmse_all'];
CORs(2:s(4)+2) = (STDs(2:s(4)+2).^2 + STDs(1)^2- RMSs(2:s(4)+2).^2)./(2.*STDs(2:s(4)+2)*STDs(1));

labels = {'EnsMean','RAW','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};

% %Removal of outliers is not applied
% q25 = prctile(RMSs(3:end),25);
% q75 = prctile(RMSs(3:end),75);
%outliers(3:s(4)+2) = RMSs(3:end)>q75+1.5*(q75-q25);
%outliers(3:s(4)+2) = RMSs(3:end)>q75*2;
% disp(['Outliers: ',labels(outliers)])
% STDs = STDs(~outliers);
% CORs = CORs(~outliers);
% RMSs = RMSs(~outliers);
% labels = labels(~outliers);

if startyear == 1980
    switch variable
        case 'tos'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.6,1.2,0.57,0.7,21) 
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.6,1.3,0.7,0.8,25)
            pretty_figure(700,650,'none','none','none','none',16);
        case 'tas'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.6,1.2,0.55,0.7,22) 
        case 'pr'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0,1.3,2.15,2.4,25)
        case 'psl'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0,1.3,1.4,1.8,19)
        case 'monmaxpr'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0,1.3,3.32,3.6,19)
        case 'monmaxtasmax'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.6,1.2,0.57,0.7,22)
        case 'monmintasmin'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.4,1.3,0.8,1,21)
        case 'zmta'
            plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end),0.9,1.2,0.4,0.45,19)
    end
else
    plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end))
end