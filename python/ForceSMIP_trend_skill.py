# %% imports
import xcdat as xc
import glob
import numpy as np
import xarray as xr
from fx import get_delta_trend, seasonal_average, global_mean

# %% Parameters
variable = 'tos'
dpath = '/p/user_pub/climate_work/pochedley1/ForceSMIP_Tier1_final/'
members = ['1B','1D','1E','1G','1J']
season = 'annual'
obs_member = '1I'
startyear = 1980
endyear = 2022
variableDict = {'monmaxtaxmax': 'tasmax',
                'monmaxpr': 'pr',
                'monmintasmin': 'tasmin',
                'zmta': 'ta'}
submission_short_names = {'RegGMST', '4th-Order-Polynomial', '10yr-Lowpass', 'LFCA', 'LFCA-2',
                          'MF-LFCA', 'MF-LFCA-2', 'LIMnMCA', 'DMDc', 'ICA-lowpass', 'LIMopt',
                          'LIMopt-filter', 'Colored-LIMnMCA', 'GPCA', 'GPCA-DA', 'RegGMST-LENSem',
                          'AllFinger', 'MLR-Forcing', 'MonthFinger', '3DUNet-Fingerprinters',
                          'SNMP-OF', 'EOF-SLR', 'LDM-SLR', 'Anchor-OPLS', 'UNet3D-LOCEAN',
                          'TrainingEM', 'RandomForest', 'EncoderDecoder', 'EnsFMP',
                          'ANN-Fingerprinters'};
seasons = xr.DataArray(data=['DJF', 'MAM', 'JJA', 'SON'],
                       dims='season',
                       coords={'season': ['DJF', 'MAM', 'JJA', 'SON']})

# %% get files
submission_files = glob.glob(dpath + 'submissions_standardized/*nc')
emean_files = glob.glob(dpath + 'ensmeans/*nc')
ref_files = glob.glob(dpath + 'Evaluation-Tier1/*nc')
if variable in variableDict.keys():
    varnam = variableDict[variable]
else:
    varnam = variable

# %% loop over members of interest
for j, member in enumerate(members):
    print(str(j+1) + ' / ' + str(len(members)) + ': ' + member)
    # get simplicity order
    if variable == 'zmta':
        simplicity_order = [np.nan, 2, 1, 12, 13, 17, 18, 14, np.nan, 10, 15, 16, 4, 8, 9, np.nan, 3, 19, np.nan, np.nan, 20, 5, 11, np.nan, 22, 21, np.nan, 6, 7, np.nan];
    else:
        simplicity_order = [26, 3, 1, 16, 17, 21, 22, 18, 8, 14, 19, 20, 7, 12, 13, 27, 5, 23, 24, 2, 28, 9, 15, 6, 30, 29, 25, 10, 11, 4];
    # get submission, ensemble mean, and reference files
    submission_file = [fn for fn in submission_files if '/' + variable + '_' + member in fn][0]
    ## NOTE: Why do we take 1A if it is obs?
    ## I assume this lets the rest of the code run and you just ignore some
    ## results for observations?
    emean_file = [fn for fn in emean_files if '/' + member + '.' + variable + '.' in fn]
    if member == obs_member:
        emean_file = [fn for fn in emean_files if '/' + '1A' + '.' + variable + '.' in fn][0]
    else:
        emean_file = emean_file[0]
    ref_file = [fn for fn in [fn for fn in ref_files if '/' + variable + '_' in fn] if '_' + member + '.' in fn][0]

    # Load Data
    dss = xc.open_dataset(submission_file)
    fields = dss['forced_component'].load()
    dse = xc.open_dataset(emean_file)
    field_emean = dse['arr_EM'].load()
    dsr = xc.open_dataset(ref_file)
    field_ref = dsr[varnam].load()

    # ensure fields are masked as appropriate
    ## Note: How do we know we should mask zero?
    fields = fields.where(fields != 0., np.nan)
    fields = fields.where(fields < 1e10, np.nan)
    field_emean = field_emean.where(field_emean != 0., np.nan)
    field_emean = field_emean.where(field_emean < 1e10, np.nan)
    field_ref = field_ref.where(field_ref != 0., np.nan)
    field_ref = field_ref.where(field_ref < 1e10, np.nan)

    # Re-order data in simplicity order
    # note: need to revisit this for zmta which has nan values
    nanmember = fields.isel(member=0).copy()
    nanmember[:] = np.nan
    tmp = []
    nm = 0
    ## Note: need to test this for zmta
    for im in simplicity_order:
        if np.isnan(im):
            nanmember.member = 'nanmember' + str(nm)
            nm += 1
            tmp.append(nanmember)
        else:
            tmp.append(fields.isel(member=im-1))
    fields = xr.concat(tmp, dim='member')

    # get monthly departures
    climref = fields.groupby('time.month').mean(dim='time')
    fields = fields.groupby('time.month') - fields.groupby('time.month').mean(dim='time')
    field_emean = field_emean.groupby('time.month') - field_emean.groupby('time.month').mean(dim='time')
    field_ref = field_ref.groupby('time.month') - field_ref.groupby('time.month').mean(dim='time')

    # create seasonal means
    if 'max' in variable:
        mode = 'max'
    elif 'min' in variable:
        mode = 'min'
    else:
        mode = 'mean'

    # Get seasonal averages
    field_ref = seasonal_average(field_ref, season, mode)
    field_emean = seasonal_average(field_emean, season, mode)
    fields_seasonal = seasonal_average(fields, season, mode)

    # get trends
    trends = get_delta_trend(fields_seasonal, startyear, endyear)
    trend_ref = get_delta_trend(field_ref, startyear, endyear)
    trend_emean = get_delta_trend(field_emean, startyear, endyear)

    # ensure common mask
    if variable in ('tos', 'zmta'):
        trends = xr.where(~np.isnan(trend_ref), trends, np.nan)
        trend_emean = xr.where(~np.isnan(trend_ref), trend_emean, np.nan)

    # pre-allocate
    if j == 0:
        rmse_ref = np.zeros(len(members))
        rmse_all = np.zeros((len(trends.member), len(members)))
        std_emean = np.zeros(len(members))
        std_ref = np.zeros(len(members))
        std_all = np.zeros((len(trends.member), len(members)))
        corr_ref = np.zeros(len(members))
        corr_all = np.zeros((len(trends.member), len(members)))

    # compute stats
    # Note: Need to deal with vertical grid in vertical averager
    rmse_ref[j] = np.sqrt(global_mean((trend_emean-trend_ref)**2))
    rmse_all[:, j] = np.sqrt(global_mean((trend_emean-trends)**2))

    std_emean[j] = np.sqrt(global_mean(trend_emean**2))
    std_ref[j] = np.sqrt(global_mean(trend_ref**2))
    std_all[:, j] = np.sqrt(global_mean((trends)**2))

    corr_ref[j] = float(global_mean(trend_emean*trend_ref/np.sqrt(global_mean(trend_emean**2)*global_mean(trend_ref**2))).values)
    corr_all[:, j] = global_mean(trend_emean*trends/np.sqrt(global_mean(trend_emean**2)*global_mean(trends**2))).values

# # %% Process Observations
# submission_file = [fn for fn in submission_files if '/' + variable + '_' + obs_member in fn][0]
# ref_file = [fn for fn in [fn for fn in ref_files if '/' + variable + '_' in fn] if '_' + obs_member + '.' in fn][0]
# ## There is some logic about choosing monmax if the variable is pr ???
# # Load data
# dss = xc.open_dataset(submission_file)
# fields = dss['forced_component'].load()
# dsr = xc.open_dataset(ref_file)
# field_ref = dsr[varnam].load()

# # ensure fields are masked as appropriate
# ## Note: How do we know we should mask zero?
# fields = fields.where(fields != 0., np.nan)
# fields = fields.where(fields < 1e10, np.nan)
# field_ref = field_ref.where(field_ref != 0., np.nan)
# field_ref = field_ref.where(field_ref < 1e10, np.nan)

# # Re-order data in simplicity order
# # note: need to revisit this for zmta which has nan values
# nanmember = fields.isel(member=0).copy()
# nanmember[:] = np.nan
# tmp = []
# nm = 0
# ## Note: need to test this for zmta
# for im in simplicity_order:
#     if np.isnan(im):
#         nanmember.member = 'nanmember' + str(nm)
#         nm += 1
#         tmp.append(nanmember)
#     else:
#         tmp.append(fields.isel(member=im-1))
# fields = xr.concat(tmp, dim='member')

# # Get seasonal averages
# field_ref = seasonal_average(field_ref, season, mode)
# fields_seasonal = seasonal_average(fields, season, mode)

# # get trends
# trends = get_delta_trend(fields_seasonal, startyear, endyear)
# trend_ref = get_delta_trend(field_ref, startyear, endyear)

# # ensure common mask
# if variable in ('tos', 'zmta'):
#     trends = xr.where(~np.isnan(trend_ref), trends, np.nan)

# %% normalized before averaging
std_ref = np.sqrt(np.mean((std_ref/std_emean)**2))
std_all = np.sqrt(np.mean((std_all/std_emean)**2, axis=1))

rmse_ref = np.sqrt(np.mean((rmse_ref/std_emean)**2))
rmse_all = np.sqrt(np.mean((rmse_all/std_emean)**2, axis=1))

std_emean = np.sqrt(np.mean(std_emean**2));

# % not used in Taylor diagram
corr_ref = np.sqrt(np.mean(corr_ref**2));
corr_all = np.sqrt(np.mean(corr_all**2, axis=1));

# %% compile statistics for tailor diagram
STDs = np.array([1, std_ref] + list(std_all))
RMSs = np.array([0, rmse_ref] + list(rmse_all))
CORs = (STDs[1:]**2 + STDs[0]**2 - RMSs[1:]**2 )/ (2*STDs[1:]*STDs[0])
CORs = np.insert(CORs, 0, 1., axis=0)
q25 = np.percentile(RMSs[2:], 25, method='midpoint')  # for some reason I get a slightly different result with np.percentile
q75 = np.percentile(RMSs[2:], 75, method='midpoint')

# %%
outliers = RMSs > q75*2  # not implemented yet
# note I changed the first label to "Correct Answer"
labels = ['Correct Answer','RAW','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']

plot_taylor(STDs, std_ref, CORs, labels, xmin=0.6, xmax=1.2, ymax=0.69, cmax=0.8, cticks=25)

