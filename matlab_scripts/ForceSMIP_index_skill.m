
% Used to make Fig. 7 of Wills et al. 2025

variable = 'tos'; % 'tos','tas','pr','psl','monmaxpr','monmaxtasmax','monmintasmin','zmta'

switch variable % Sahel precip uses MJJAS, Aleutian Low SLP uses DJF, otherwise annual mean
    case 'pr'
        season = [5 6 7 8 9]; 
    case 'psl'
        season = [12 1 2];
    otherwise
        season = 1:12;
end

submission_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/submissions-Tier1-standardized-estimates');
emean_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/ensmeans-Tier1');
ref_files = get_files('/Users/rjnglin/Data/ForceSMIP/ForceSMIP_Tier1_final/Evaluation-Tier1');

member = {'1B','1D','1E','1G','1J'}; % unseen models '1B','1D','1E','1G','1J' in ForceSMIP Tier 1
obs = {'1I'};

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

rmse_ref = zeros(1,length(member)); std_emean = rmse_ref; std_ref = rmse_ref; corr_ref = rmse_ref;
rmse_all = zeros(30,length(member)); std_all = rmse_all; corr_all = rmse_all;

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
    if strcmp(variable,'pr')
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
    years = floor(1:1/12:nyr+0.99);

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

    % compute indices
    switch variable
        case 'tas' % GMST
            index_emean = global_mean(lon,lat,field_emean);
            index_ref = global_mean(lon,lat,field_ref);
            index_all = global_mean(lon,lat,fields_seasonal);
        case 'monmaxtasmax' % Europe land
            land = mask; land(~isnan(mask)) = nan; land(isnan(mask)) = 1;
            index_emean = mean_in_a_box(lon,lat,land.*field_emean,[40 55],[0 40]);
            index_ref = mean_in_a_box(lon,lat,land.*field_ref,[40 55],[0 40]);
            index_all = mean_in_a_box(lon,lat,land.*fields_seasonal,[40 55],[0 40]);
        case 'monmaxpr' % Boulder
%             index_emean = global_mean(lon,lat,field_emean);
%             index_ref = global_mean(lon,lat,field_ref);
%             index_all = global_mean(lon,lat,fields_seasonal);
            index_emean = rmean(squeeze(field_emean(closest(lon,360-105.27),closest(lat,40.01),:))',30);
            index_ref = rmean(squeeze(field_ref(closest(lon,360-105.27),closest(lat,40.01),:))',30);
            for i = 1:size(fields,4)
                index_all(:,i) = rmean(squeeze(fields_seasonal(closest(lon,360-105.27),closest(lat,40.01),:,i)),30);
            end
        case 'tos' % Nino34 minus GMSST
%             index_emean = rmean(mean_in_a_box(lon,lat,field_emean,[-5 5],[170 240])-global_mean(lon,lat,field_emean),10);
%             index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[-5 5],[170 240])-global_mean(lon,lat,field_ref),10);
%             for i = 1:size(fields,4)
%                 index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[-5 5],[170 240])-global_mean(lon,lat,fields_seasonal(:,:,:,i)),10);
%             end
            % NASSTI minus GMSST
            index_emean = rmean(mean_in_a_box(lon,lat,field_emean,[0 60],[280 360])-global_mean(lon,lat,field_emean),10);
            index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[0 60],[280 360])-global_mean(lon,lat,field_ref),10);
            for i = 1:size(fields,4)
                index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[0 60],[280 360])-global_mean(lon,lat,fields_seasonal(:,:,:,i)),10);
            end
        case 'pr' % Sahel rainfall
            index_emean = rmean(mean_in_a_box(lon,lat,field_emean,[10 20],[340 10]),10);
            index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[10 20],[340 10]),10);
            for i = 1:size(fields,4)
                index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[10 20],[340 10]),10);
            end
        case 'psl' % Aleutian low (NPI)
            index_emean = rmean(mean_in_a_box(lon,lat,field_emean,[30 65],[160 220]),10);
            index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[30 65],[160 220]),10);
            for i = 1:size(fields,4)
                index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[30 65],[160 220]),10);
            end
            % NAO (with units)
            %index_emean = squeeze(field_emean(closest(lon,360-21.8174),closest(lat,64.1265),:)-field_emean(closest(lon,360-9.1393),closest(lat,38.7223),:))';
            %index_ref = squeeze(field_ref(closest(lon,360-21.8174),closest(lat,64.1265),:)-field_ref(closest(lon,360-9.1393),closest(lat,38.7223),:))';
            %index_all = squeeze(fields_seasonal(closest(lon,360-21.8174),closest(lat,64.1265),:,:)-fields_seasonal(closest(lon,360-9.1393),closest(lat,38.7223),:,:));
    end

    index_emean = index_emean-mean(index_emean);
    index_ref = index_ref-mean(index_ref);
    index_all = index_all-mean(index_all,1);

    % compute skll scores
    for i = 1:s(4)
        if ~strcmp(variable,'zmta')
            rmse_ref(j) = sqrt(mean((index_emean-index_ref).^2));
            rmse_all(i,j) = sqrt(mean((index_emean-index_all(:,i)').^2));
            std_emean(j) = sqrt(mean(index_emean.^2));
            std_ref(j) = sqrt(mean(index_ref.^2));
            std_all(i,j) = sqrt(mean(index_all(:,i).^2));
            corr_ref(j) = mean(index_emean.*index_ref)./(std_emean(j).*std_ref(j));
            corr_all(i,j) = mean(index_emean.*index_all(:,i)')./(std_emean(j).*std_all(i,j));
        else
            warning('No zmta index has been written into this code.')
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

% compute indices
switch variable
    case 'tas' % GMST
        index_ref = global_mean(lon,lat,field_ref);
        index_all = global_mean(lon,lat,fields_seasonal);
    case 'monmaxtasmax' % Europe land
        land = mask; land(~isnan(mask)) = nan; land(isnan(mask)) = 1;
        index_emean = mean_in_a_box(lon,lat,land.*field_emean,[40 55],[0 40]);
        index_ref = mean_in_a_box(lon,lat,land.*field_ref,[40 55],[0 40]);
        index_all = mean_in_a_box(lon,lat,land.*fields_seasonal,[40 55],[0 40]);
    case 'monmaxpr' % Boulder
%         index_ref = global_mean(lon,lat,field_ref);
%         index_all = global_mean(lon,lat,fields_seasonal);
        index_ref = rmean(squeeze(field_ref(closest(lon,360-105.27),closest(lat,40.01),:))',30);
        for i = 1:size(fields,4)
            index_all(:,i) = rmean(squeeze(fields_seasonal(closest(lon,360-105.27),closest(lat,40.01),:,i)),30);
        end
    case 'tos' % Nino34 minus GMSST
%         index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[-5 5],[170 240])-global_mean(lon,lat,field_ref),10);
%         for i = 1:size(fields,4)
%             index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[-5 5],[170 240])-global_mean(lon,lat,fields_seasonal(:,:,:,i)),10);
%         end
        index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[0 60],[280 360])-global_mean(lon,lat,field_ref),10);
        for i = 1:size(fields,4)
            index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[0 60],[280 360])-global_mean(lon,lat,fields_seasonal(:,:,:,i)),10);
        end
    case 'pr' % Sahel rainfall
        index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[10 20],[340 10]),10);
        for i = 1:size(fields,4)
            index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[10 20],[340 10]),10);
        end
    case 'psl' % Aluetian Low (NPI)
        index_ref = rmean(mean_in_a_box(lon,lat,field_ref,[30 65],[160 220]),10);
        for i = 1:size(fields,4)
            index_all(:,i)  = rmean(mean_in_a_box(lon,lat,fields_seasonal(:,:,:,i),[30 65],[160 220]),10);
        end
end

clear fields_seasonal

%% averaging skill metrics

index_ref = index_ref-mean(index_ref);
index_all = index_all-mean(index_all,1);

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
q25 = prctile(RMSs(3:end),25);
q75 = prctile(RMSs(3:end),75);
%outliers(3:s(4)+2) = RMSs(3:end)>q75+1.5*(q75-q25);
outliers(3:s(4)+2) = RMSs(3:end)>q75*2;
labels = {'EnsMean','RAW','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
% disp(['Outliers: ',labels(outliers)])
% STDs = STDs(~outliers);
% CORs = CORs(~outliers);
% RMSs = RMSs(~outliers);
% labels = labels(~outliers);
plot_Taylor(STDs(2:end),1,CORs(2:end),labels(2:end))