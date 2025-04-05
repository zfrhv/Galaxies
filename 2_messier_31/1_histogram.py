import plotly.graph_objects as go
import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack

# import csv results
results = Table.read('andromeda_square_data.csv', format='csv')

# filter results
results = results[results['ra_error'] > 0.1] # idk why this Gaia is so uncertain but smaller than 0.1 looks like noise
results = results[results['dec_error'] > 0.1]
results = results[~np.isnan(results['bp_rp'])]
results = results[~results['bp_rp'].mask]
results = results[~np.isnan(results['bp_g'])]
results = results[~results['bp_g'].mask]
results = results[~np.isnan(results['g_rp'])]
results = results[~results['g_rp'].mask]

# sort stars for histograme
results.sort('ra')

screen_parts = 4
histogram_parts = 1000

histogram_height = histogram_parts * screen_parts * 1/10000

# split X axis (ra) into screen_parts
min_ra = np.min(results['ra'])
max_ra = np.max(results['ra']) + 0.000001
bins_x = np.linspace(min_ra, max_ra, screen_parts+1)

split_results = []
for i in range(screen_parts):
    lower, upper = bins_x[i], bins_x[i+1]
    subset = results[(results['ra'] >= lower) & (results['ra'] < upper)]
    split_results.append(subset)

# split Y axis (dec) into screen_parts
min_dec = np.min(results['dec'])
max_dec = np.max(results['dec']) + 0.000001
bins_y = np.linspace(min_dec, max_dec, screen_parts+1)

for i in range(len(split_results)):
    sub_results = split_results[i]
    split_sub_results = []
    for j in range(screen_parts):
        lower, upper = bins_y[j], bins_y[j+1]
        subset = sub_results[(sub_results['dec'] >= lower) & (sub_results['dec'] < upper)]
        split_sub_results.append(subset)
    split_results[i] = split_sub_results

# create histogram inside each bin
noise_histogram_bp_rp = None
noise_histogram_bp_g = None
noise_histogram_g_rp = None
min_bp_rp = np.min(results['bp_rp'])
max_bp_rp = np.max(results['bp_rp']) + 0.000001
min_bp_g = np.min(results['bp_g'])
max_bp_g = np.max(results['bp_g']) + 0.000001
min_g_rp = np.min(results['g_rp'])
max_g_rp = np.max(results['g_rp']) + 0.000001
bp_rp_x = []
bp_rp_y = []
bp_g_x = []
bp_g_y = []
g_rp_x = []
g_rp_y = []
for i in range(screen_parts):
    for j in range(screen_parts):
        table = split_results[i][j]
        start_x = bins_x[i]
        end_x = bins_x[i+1]
        start_y = bins_y[j]
        end_y = bins_y[j+1]

        histogram_g = np.zeros(histogram_parts+1, dtype=float)
        histogram_bp = np.zeros(histogram_parts+1, dtype=float)
        histogram_rp = np.zeros(histogram_parts+1, dtype=float)
        section_size = (max_ra - min_ra) / screen_parts / histogram_parts
        for result in table:
            bp_rp = result['bp_rp']
            bp_rp_ratio = (bp_rp - min_bp_rp) / (max_bp_rp - min_bp_rp)
            section = math.floor(bp_rp_ratio * histogram_parts)
            offset_x = section * section_size
            offset_y = histogram_g[section] * (section_size * histogram_height)
            histogram_g[section] += 1
            bp_rp_x.append(start_x + offset_x)
            bp_rp_y.append(start_y + offset_y)

            bp_g = result['bp_g']
            bp_g_ratio = (bp_g - min_bp_g) / (max_bp_g - min_bp_g)
            section = math.floor(bp_g_ratio * histogram_parts)
            offset_x = section * section_size + 1/3*section_size
            offset_y = histogram_bp[section] * (section_size * histogram_height)
            histogram_bp[section] += 1
            bp_g_x.append(start_x + offset_x)
            bp_g_y.append(start_y + offset_y)

            g_rp = result['g_rp']
            g_rp_ratio = (g_rp - min_g_rp) / (max_g_rp - min_g_rp)
            section = math.floor(g_rp_ratio * histogram_parts)
            offset_x = section * section_size + 2/3*section_size
            offset_y = histogram_rp[section] * (section_size * histogram_height)
            histogram_rp[section] += 1
            g_rp_x.append(start_x + offset_x)
            g_rp_y.append(start_y + offset_y)

        # use the right bottom histogram for noise inspection
        if screen_parts == 3 and i == screen_parts-1 and j == 0:
            noise_histogram_bp_rp = histogram_g
            noise_histogram_bp_g = histogram_bp
            noise_histogram_g_rp = histogram_rp
        # use the left up histogram for noise inspection
        if screen_parts == 4 and i == 0 and j == screen_parts-1:
            noise_histogram_bp_rp = histogram_g
            noise_histogram_bp_g = histogram_bp
            noise_histogram_g_rp = histogram_rp



# merge the results back into one table
flattened_results = []
for sublist in split_results:
    for table in sublist:
        flattened_results.append(table)
split_results = flattened_results

results = vstack(split_results)

# # get useful data
# ra = np.array(results['ra'])
# dec = np.array(results['dec'])
# bp_rp_mean_mag = np.array(results['bp_rp_mean_mag'])
# bp_g_mean_mag = np.array(results['bp_g_mean_mag'])
# g_rp_mean_mag = np.array(results['g_rp_mean_mag'])
# bp_rp_mean_mag_error = np.array(results['bp_rp_mean_mag_error'])
# bp_g_mean_mag_error = np.array(results['bp_g_mean_mag_error'])
# g_rp_mean_mag_error = np.array(results['g_rp_mean_mag_error'])

# Create a 2D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scattergl(
    x=bp_rp_x,
    y=bp_rp_y,
    mode='markers',
    marker=dict(
        size=1,
        color='green',
        opacity=1
    ),
))

fig.add_trace(go.Scattergl(
    x=bp_g_x,
    y=bp_g_y,
    mode='markers',
    marker=dict(
        size=1,
        color='lightblue',
        opacity=1
    ),
))

fig.add_trace(go.Scattergl(
    x=g_rp_x,
    y=g_rp_y,
    mode='markers',
    marker=dict(
        size=1,
        color='red',
        opacity=1
    ),
))

# Labels and title
fig.update_layout(
    title="2D Plot of Star Positions (RA vs DEC)",
    xaxis_title='Right Ascension (RA) [degrees]',
    yaxis_title='Declination (DEC) [degrees]',
    plot_bgcolor='black',
    paper_bgcolor='black',
    font=dict(color='white'),
    xaxis=dict(
        showgrid=False,
        zeroline=False,
        scaleanchor='y',
    ),
    yaxis=dict(
        showgrid=False,
        zeroline=False,
    ),
)

# Show the plot
fig.show()


# used bottom right square size of 1/3 out of the screen (screen_parts = 3)
# it doesnt has any special stuff and only noise stars
if screen_parts == 3 and histogram_parts == 1000:
    area = ((max_ra - min_ra)/screen_parts * (max_dec - min_dec)/screen_parts)
    np.save('noise_histograms/noise_histogram_bp_rp_right_down.npy', noise_histogram_bp_rp)
    np.save('noise_histograms/noise_histogram_bp_g_right_down.npy', noise_histogram_bp_g)
    np.save('noise_histograms/noise_histogram_g_rp_right_down.npy', noise_histogram_g_rp)
    np.save('noise_histograms/noise_area_right_down.npy', area)
if screen_parts == 4 and histogram_parts == 1000:
    area = ((max_ra - min_ra)/screen_parts * (max_dec - min_dec)/screen_parts)
    np.save('noise_histograms/noise_histogram_bp_rp_left_up.npy', noise_histogram_bp_rp)
    np.save('noise_histograms/noise_histogram_bp_g_left_up.npy', noise_histogram_bp_g)
    np.save('noise_histograms/noise_histogram_g_rp_left_up.npy', noise_histogram_g_rp)
    np.save('noise_histograms/noise_area_left_up.npy', area)