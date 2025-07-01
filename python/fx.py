# imports
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# helper functions
def get_delta_trend(data, startyear, endyear):
    data = data.copy()
    data_period = data.sel(year=slice(startyear, endyear))
    data_trend = data_period.polyfit(dim='year', deg=1).polyfit_coefficients.isel(degree=0)
    tdelta = data_period.year.values[-1] - data_period.year.values[0]
    data_delta = data_trend * tdelta
    return data_delta

def seasonal_average(data, season, mode):
    data = data.copy()
    # if season is annual, do calculation and return result
    if season == 'annual':
        return data.groupby(data.time.dt.year).mean("time")
    # get vectors of years / months
    months = data.time.dt.month
    years = data.time.dt.year
    # get original time vector (to be modified if needed for DJF)
    time = data.time.copy()
    # set December values to year+1 (so that they can be averaged
    # with JF of the following year)
    if season == 'DJF':
        # make sure time series starts with January
        # we take JF as the first DJF average
        if time.values[0].month != 1:
            ValueError("Expecting array to start in January")
        # find December months
        I = np.where(months == 12)[0]
        # December time axis
        dectime = time[I]
        # Add one year to each December time value
        ctype = type(dectime.values[0])
        newdectime = [ctype(t.year+1, t.month, t.day, t.hour, t.minute, t.second) for t in dectime.values]
        # Update time with modified Decembers
        time.values[I] = newdectime
        # Update time axis in original dataarray for grouping calculations
        data['time'] = time
    # get seasonal groupings
    data_season = data.sel(time=data.time.dt.season==season)
    # take min/mean/max by year
    if mode == 'mean':
        data_out = data_season.groupby(data_season.time.dt.year).mean("time")
    elif mode == 'min':
        data_out = data_season.groupby(data_season.time.dt.year).min("time")
    elif mode == 'max':
        data_out = data_season.groupby(data_season.time.dt.year).min("time")
    # remove extra year from DJF calculation
    data_out = data_out.sel(year=list(set(years.values)))
    return data_out

def global_mean(data):
    data = data.copy()
    lat = data.lat.values
    diffs = lat[1:] - lat[0:-1]
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])
    lower_bounds = np.sin(np.radians(lat - diffs[:-1] * 0.5))
    upper_bounds = np.sin(np.radians(lat + diffs[1:] * 0.5))
    weights = np.abs(upper_bounds - lower_bounds)
    weights = xr.DataArray(weights,
                           dims=['lat'],
                           coords={'lat': data.lat})
    daw = data.weighted(weights)
    return daw.mean(dim=('lat', 'lon'))

def plot_taylor(STDs, STD_ref, CORs, labels, xmin=None, xmax=None, ymax=None, cmax=None, cticks=None):
    plt.figure(figsize=(6, 6.5))
    alpha = np.arccos(CORs)
    x = STDs * np.cos(alpha)
    y = STDs * np.sin(alpha)
    labels = np.array(labels)[~np.isnan(x)]
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if xmin is None:
        ymax = np.ceil(np.max(y*10)/10)
        xmax = np.ceil(np.max(x*10)/10)
        xmin = np.min([1, np.floor(np.min(x*10)/10)])
    xi = np.arange(xmin, xmax+0.01, 0.01)
    yi = np.arange(0, ymax+0.01, 0.01)
    [Y,X] = np.meshgrid(yi,xi)
    ci = np.sqrt((X-1)**2+Y**2)
    ci2 = np.sqrt(X**2+Y**2)
    if cmax is None:
        cmax = np.ceil(np.max(np.max(ci))*10)/10;
        cticks = 21
    if labels is not None:
        if xmax < 2:
            dx = xmax/200
            dy = ymax/200
        else:
            dx = xmax/100
            dy = ymax/100
    # %%
    im = plt.contourf(xi, yi, ci.T, np.linspace(0, cmax, cticks), cmap=plt.cm.YlOrRd)
    plt.plot(x[1:], y[1:], 'kd', markersize=5)
    for i in range(len(x)):
        if ((y[i] < ymax) & (x[i] > xmin) & (x[i] < xmax)):
            plt.text(x[i]+dx, y[i]+dy, labels[i], fontsize=8)
    plt.plot(x[0], y[0], marker='p', color='k', markersize=5)
    rvalues = [-0.99, -0.95] + list(np.arange(-0.9, 0, 0.1)) + list(np.arange(0.1, 0.91, 0.1)) + [0.95, 0.99]
    yticks = []
    yticklabels = []
    for r in rvalues:
        ycorr = xi * np.tan(np.arccos(r))
        plt.plot(xi, ycorr,'k')
        yticks.append(np.interp(xmax, xi, ycorr))
        yticklabels.append(r)
    ax = plt.gca()
    plt.yticks(yticks, labels=yticklabels)
    plt.plot([0, 0],[0, ymax], 'k', linewidth=1.5)
    ax.yaxis.tick_right() # Moves ticks and labels to the right
    ax.yaxis.set_label_position("right") # Optional: move the y-axis label to the right as well
    if xmax < 3:
        plt.contour(xi, yi, ci2.T, np.arange(0.2, 10.1, 0.2), colors=['k'])
    else:
        plt.contour(xi, yi, ci2.T, np.arange(0.5, 10.1, 0.5), colors=['k'])
    plt.contour(xi, yi, ci2.T, [1], colors=['k'], linewidth=1.5)
    # add white line
    xi_hi = np.arange(xmin, xmax, 0.001)
    yi_hi = np.arange(0, ymax, 0.001)
    xi_him = np.tile(np.expand_dims(xi_hi, axis=1), (1, len(yi_hi)))
    yi_him = np.tile(np.expand_dims(yi_hi, axis=1).T, (len(xi_hi), 1))
    fred_rmse = 1 - np.sqrt((xi_him - 1)**2 + yi_him**2 ) / np.sqrt((x[1]-1)**2 + y[1]**2)
    fred_pcorr = 1 - (xi_him / np.sqrt(xi_him**2 + yi_him**2)) / (x[1] / np.sqrt(x[1]**2 + y[1]**2))
    plt.contour(xi_hi, yi_hi, (fred_rmse>fred_pcorr).T, linestyles=':', colors='white', linewidth=1.5)
    # hold on; contour(xi_hi',yi_hi,(fred_rmse>fred_pcorr)',[0.5 0.5],':','color','w','linewidth',1.5);
    cb = plt.colorbar(im, orientation='horizontal', pad=0.05)
    cb.set_label('RMSE (Normalized)')
    plt.xlim(xmin, xmax)
    plt.ylim(0., ymax)
    plt.tight_layout()
    plt.show()
