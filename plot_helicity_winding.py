import scipy.io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import requests

def read_data(file_dir, file_name, start_time):
    
    data = pd.read_csv(file_dir + file_name)
    print('data head: ', data.head(0))

    times= data['times(seconds)']
    vars = data.iloc[:, 1:].values
    
    date_times = [start_time + timedelta(seconds=float(seconds)) for seconds in times]
    
    return date_times, vars

def get_cme_source_location(location, lat_min, lat_max, lon_min, lon_max):
    
    if not location:
        return None, None, False
    # Converting signed latitude
    if 'N' in location:
        lat_index = location.index('N')
        lat = float(location[lat_index + 1:lat_index + 3])
    elif 'S' in location:
        lat_index = location.index('S')
        lat = -float(location[lat_index + 1:lat_index + 3])
    # Converting signed longitude 
    if 'E' in location:
        lon_index = location.index('E')
        lon = float(location[lon_index + 1:lon_index + 3])
    elif 'W' in location:
        lon_index = location.index('W')
        lon = -float(location[lon_index + 1:lon_index + 3])
    # Judging whether CME orginates from the active region
    in_AR = False
    if lat is not None and lon is not None:
        in_AR = lat_min <= lat <= lat_max and lon_min <= lon <= lon_max
        
    return lat, lon, in_AR

def get_cme_data(startDate, endDate, noaa_number):
    
    base_url = 'https://api.nasa.gov/DONKI/CME'
    params = {
        'startDate': startDate,
        'endDate': endDate,
        'api_key': api_key
    }
    if noaa_number == 13664:
        lat_min, lat_max, lon_min, lon_max = -30, -10, -35, 5
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        filtered_cme = []
        for cme in data:
            cme_type = cme.get('cmeAnalyses', [{}])[0].get('type')
            # Judging CME source location
            cme_location = cme.get('sourceLocation')
            lat, lon, in_AR = get_cme_source_location(cme_location, lat_min, lat_max, lon_min, lon_max)
            # Recording CME info
            if cme_type in ['C', 'O', 'R', 'ER'] and in_AR:
                filtered_cme.append({
                    'startTime': cme.get('startTime'),
                    'type': cme_type
                })
        return filtered_cme
    except requests.RequestException as e:
        print(f'Request error: {e}')
    except ValueError:
        print('Failed to parse JSON data')
    return []

def get_flare_data(startDate, endDate, noaa_number):
    base_url = 'https://api.nasa.gov/DONKI/FLR'
    params = {
        'startDate': startDate,
        'endDate': endDate,
        'api_key': api_key
    }
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        filtered_flare = []
        for flare in data:
            flare_class = flare.get('classType')
            if flare_class and (flare_class.startswith('M') or flare_class.startswith('X')) \
                and flare.get('activeRegionNum') == 13664:
                filtered_flare.append({
                    'beginTime': flare.get('beginTime'),
                    'class': flare_class
                })
        return filtered_flare
    except requests.RequestException as e:
        print(f'Request error: {e}')
    except ValueError:
        print('Failed to parse JSON data')
    return []

def plot_figure(date_times, wind, hel, windBraid, helBraid, \
    wind_int, hel_int, windBraid_int, helBraid_int, \
    cme_time_lst, cme_type_lst, flare_time_lst, flare_type_lst):
    
    cme_color_map = {
        'C': 'orange',
        'O': 'red',
        'R': 'black',
        'ER': 'black'
    }
    flare_color_map = {
        'M': 'blue',
        'X': 'purple'
    }

    def plot_event_lines(ax, event_type, time_lst, type_lst, color_map, linestyle='--'):
        
        for event_time, event_type in zip(time_lst, type_lst):
            if event_time and event_type:
                try:
                    event_datetime = datetime.strptime(event_time, '%Y-%m-%dT%H:%MZ')
                    color = color_map.get(event_type[0], 'gray')
                    ax.axvline(x=event_datetime, color=color, linewidth=2, linestyle=linestyle, alpha=0.7)
                except ValueError:
                    print(f'Failed to parse time: {event_time}')
            else:
                print('There is no event!')
    
    # Counting event numbers
    c_cme_count = cme_type_lst.count('C')
    o_cme_count = cme_type_lst.count('O')
    r_er_cme_count = cme_type_lst.count('R') + cme_type_lst.count('ER')
    m_flare_count = sum(1 for t in flare_type_lst if t and t.startswith('M'))
    x_flare_count = sum(1 for t in flare_type_lst if t and t.startswith('X'))
    
    plt.figure(figsize=(8,7))
    plt.rcParams.update({'font.size': 14})
    xmin = np.min(date_times)
    xmax = np.max(date_times)

    # ax1 = plt.subplot(4,1,1)
    # plt.plot(date_times, helBraid, label='braid')
    # plt.plot(date_times, hel-helBraid, label='emerge')
    # plt.plot(date_times, hel, label='total')
    # # plt.xticks([])
    # plt.xlim([xmin, xmax])
    # plt.ylabel(r'$\frac{dH}{dt} [G^2 km^4/s]$')
    # # plt.yscale('symlog')
    # plot_event_lines(ax1, 'cme', cme_time_lst, cme_type_lst, cme_color_map)
    # plot_event_lines(ax1, 'flare', flare_time_lst, flare_type_lst, flare_color_map)
    # date_format = mdates.DateFormatter('%m-%d-%H')
    # ax1.xaxis.set_major_formatter(date_format)

    # ax2 = plt.subplot(4,1,2, sharex=ax1)
    # plt.plot(date_times, windBraid, label='braid')
    # plt.plot(date_times, wind-windBraid, label='emerge')
    # plt.plot(date_times, wind, label='total')
    # # plt.xticks([])
    # plt.xlim([xmin, xmax])
    # plt.ylabel(r'$\frac{dL}{dt} [km^4/s]$')
    # # plt.yscale('symlog')
    # plot_event_lines(ax2, 'cme', cme_time_lst, cme_type_lst, cme_color_map)
    # plot_event_lines(ax2, 'flare', flare_time_lst, flare_type_lst, flare_color_map)

    # ax3 = plt.subplot(4,1,3, sharex=ax1)
    # plt.plot(date_times, helBraid_int, label='braid')
    # plt.plot(date_times, hel_int-helBraid_int, label='emerge')
    # plt.plot(date_times, hel_int, label='total')
    # # plt.xticks([])
    # plt.xlim([xmin, xmax])
    # plt.ylabel(r'$H [G^2 km^4]$')
    # plt.yscale('symlog')
    # plot_event_lines(ax3, 'cme', cme_time_lst, cme_type_lst, cme_color_map)
    # plot_event_lines(ax3, 'flare', flare_time_lst, flare_type_lst, flare_color_map)

    # ax4 = plt.subplot(4,1,4, sharex=ax1)
    ax4 = plt.subplot(1,1,1)
    plt.hlines(y=0, xmin=xmin, xmax=xmax, colors='gray')
    plt.plot(date_times, windBraid_int, color='blue', label=r'$L_{braid}$')
    plt.plot(date_times, wind_int-windBraid_int, color='green', label=r'$L_{emerge}$')
    plt.plot(date_times, wind_int, color='black', label=r'$L$')
    # plt.legend(prop={'size': 12})
    plt.xlim([xmin, xmax])
    plt.xlabel('Date')
    plt.ylabel(r'Winding Accumulation [km$^4$]')
    # plt.yscale('symlog')
    plot_event_lines(ax4, 'cme', cme_time_lst, cme_type_lst, cme_color_map, linestyle = '--')
    plot_event_lines(ax4, 'flare', flare_time_lst, flare_type_lst, flare_color_map, linestyle = ':')
    date_format = mdates.DateFormatter('%m-%d')
    plt.gca().xaxis.set_major_formatter(date_format)
    
    # Add legend for events
    legend_elements = [
        Line2D([0], [0], color='blue', lw=2, label=r'$L_{braid}$'),
        Line2D([0], [0], color='green', lw=2, label=r'$L_{emerge}$'),
        Line2D([0], [0], color='black', lw=2, label=r'$L$'),
        Line2D([0], [0], color='orange', lw=2, label=f'C-class CME ({c_cme_count})', linestyle='--'),
        Line2D([0], [0], color='red', lw=2, label=f'O-class CME ({o_cme_count})', linestyle='--'),
        # Line2D([0], [0], color='black', lw=2, label=f'R/ER-class CME ({r_er_cme_count})', linestyle='--'),
        # Line2D([0], [0], color='blue', lw=2, label=f'M-class Flare ({m_flare_count})', linestyle='--'),
        Line2D([0], [0], color='purple', lw=2, label=f'X-class Flare ({x_flare_count})', linestyle=':')
    ]
    # ax1.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.3), ncol=5, prop={'size': 12})
    ax4.legend(handles=legend_elements, loc='upper left', ncol=1, prop={'size': 12})

    # plt.subplots_adjust(hspace=0.05)

    plt.show()
    return

## Reading data
file_dir = 'E:/Research/Work/202405_solar_storm/helicity/run_3/'
rate_file = 'raw_series.csv'
accum_file = 'accum.csv'
start_time = datetime(2024, 5, 4, 0, 0, 0)
date_times, rate_vars = read_data(file_dir, rate_file, start_time)
date_times, accum_vars = read_data(file_dir, accum_file, start_time)

## Extracting data
windCur, helCur = rate_vars[:,0], rate_vars[:,1]
windPot, helPot = rate_vars[:,2], rate_vars[:,3]
windBraid, helBraid = rate_vars[:,4], rate_vars[:,5]
wind, hel = rate_vars[:,6], rate_vars[:,7]
dL, dH = rate_vars[:,8], rate_vars[:,9]

windCur_int, helCur_int = accum_vars[:,0], accum_vars[:,1]
windPot_int, helPot_int = accum_vars[:,2], accum_vars[:,3]
windBraid_int, helBraid_int = accum_vars[:,4], accum_vars[:,5]
wind_int, hel_int = accum_vars[:,6], accum_vars[:,7]
dL_int, dH_int = accum_vars[:,8], accum_vars[:,9]

## Getting DONKI list
api_key = 'ogX37BONEfKdzMUIdgO17Nor0rch6vKL4eH0WU1A'
startDate = '2024-05-04'
endDate = '2024-05-09'
noaa_number = 13664
cme_data = get_cme_data(startDate, endDate, noaa_number)
flare_data = get_flare_data(startDate, endDate, noaa_number)

## Extracting eruption events
cme_time_lst = []
cme_type_lst = []
for cme in cme_data:
    cme_time = cme.get('startTime')
    cme_type = cme.get('type')
    cme_time_lst.append(cme_time)
    cme_type_lst.append(cme_type)

flare_time_lst = []
flare_type_lst = []
for flare in flare_data:
    flare_time = flare.get('beginTime')
    flare_type = flare.get('class')
    if flare_type[0] == 'X':
        flare_time_lst.append(flare_time)
        flare_type_lst.append(flare_type)

## Plotting figure
plot_figure(date_times, wind, hel, windBraid, helBraid, \
    wind_int, hel_int, windBraid_int, helBraid_int, \
    cme_time_lst, cme_type_lst, flare_time_lst, flare_type_lst)

db