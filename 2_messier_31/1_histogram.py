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
results = results[~np.isnan(results['phot_g_mean_mag'])]
results = results[~results['phot_g_mean_mag'].mask]
results = results[~np.isnan(results['phot_bp_mean_mag'])]
results = results[~results['phot_bp_mean_mag'].mask]
results = results[~np.isnan(results['phot_rp_mean_mag'])]
results = results[~results['phot_rp_mean_mag'].mask]

# sort stars for histograme
results.sort('ra')

screen_parts = 4
histogram_parts = 100

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
noise_histogram_g = None
noise_histogram_bp = None
noise_histogram_rp = None
min_phot_g = np.min(results['phot_g_mean_mag'])
max_phot_g = np.max(results['phot_g_mean_mag']) + 0.000001
min_phot_bp = np.min(results['phot_bp_mean_mag'])
max_phot_bp = np.max(results['phot_bp_mean_mag']) + 0.000001
min_phot_rp = np.min(results['phot_rp_mean_mag'])
max_phot_rp = np.max(results['phot_rp_mean_mag']) + 0.000001
g_x = []
g_y = []
bp_x = []
bp_y = []
rp_x = []
rp_y = []
for i in range(screen_parts):
    for j in range(screen_parts):
        table = split_results[i][j]
        start_x = bins_x[i]
        end_x = bins_x[i+1]
        start_y = bins_y[j]
        end_y = bins_y[j+1]

        histogram_g = np.zeros(histogram_parts+1)
        histogram_bp = np.zeros(histogram_parts+1)
        histogram_rp = np.zeros(histogram_parts+1)
        section_size = (max_ra - min_ra) / screen_parts / histogram_parts
        for result in table:
            phot_g = result['phot_g_mean_mag']
            phot_g_ratio = (phot_g - min_phot_g) / (max_phot_g - min_phot_g)
            section_g = math.floor(phot_g_ratio * histogram_parts)
            offset_x_g = section_g * section_size
            offset_y_g = histogram_g[section_g] * (section_size * histogram_height)
            histogram_g[section_g] += 1
            g_x.append(start_x + offset_x_g)
            g_y.append(start_y + offset_y_g)

            phot_bp = result['phot_bp_mean_mag']
            phot_bp_ratio = (phot_bp - min_phot_bp) / (max_phot_bp - min_phot_bp)
            section_bp = math.floor(phot_bp_ratio * histogram_parts)
            offset_x_bp = section_bp * section_size + 1/3*section_size
            offset_y_bp = histogram_bp[section_bp] * (section_size * histogram_height)
            histogram_bp[section_bp] += 1
            bp_x.append(start_x + offset_x_bp)
            bp_y.append(start_y + offset_y_bp)

            phot_rp = result['phot_rp_mean_mag']
            phot_rp_ratio = (phot_rp - min_phot_rp) / (max_phot_rp - min_phot_rp)
            section_rp = math.floor(phot_rp_ratio * histogram_parts)
            offset_x_rp = section_rp * section_size + 2/3*section_size
            offset_y_rp = histogram_rp[section_rp] * (section_size * histogram_height)
            histogram_rp[section_rp] += 1
            rp_x.append(start_x + offset_x_rp)
            rp_y.append(start_y + offset_y_rp)

        # use the right bottom histogram for noise inspection
        if screen_parts == 3 and i == screen_parts-1 and j == 0:
            noise_histogram_g = histogram_g
            noise_histogram_bp = histogram_bp
            noise_histogram_rp = histogram_rp
        # use the left up histogram for noise inspection
        if screen_parts == 4 and i == 0 and j == screen_parts-1:
            noise_histogram_g = histogram_g
            noise_histogram_bp = histogram_bp
            noise_histogram_rp = histogram_rp



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
# phot_g_mean_mag = np.array(results['phot_g_mean_mag'])
# phot_bp_mean_mag = np.array(results['phot_bp_mean_mag'])
# phot_rp_mean_mag = np.array(results['phot_rp_mean_mag'])
# phot_g_mean_mag_error = np.array(results['phot_g_mean_mag_error'])
# phot_bp_mean_mag_error = np.array(results['phot_bp_mean_mag_error'])
# phot_rp_mean_mag_error = np.array(results['phot_rp_mean_mag_error'])

# Create a 2D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scattergl(
    x=g_x,
    y=g_y,
    mode='markers',
    marker=dict(
        size=1,
        color='green',
        opacity=1
    ),
))

fig.add_trace(go.Scattergl(
    x=bp_x,
    y=bp_y,
    mode='markers',
    marker=dict(
        size=1,
        color='lightblue',
        opacity=1
    ),
))

fig.add_trace(go.Scattergl(
    x=rp_x,
    y=rp_y,
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
if screen_parts == 3 and histogram_parts == 100:
    area = ((max_ra - min_ra)/screen_parts * (max_dec - min_dec)/screen_parts)
    noise_histogram_g = np.array(noise_histogram_g) / area
    noise_histogram_bp = np.array(noise_histogram_bp) / area
    noise_histogram_rp = np.array(noise_histogram_rp) / area
    np.save('noise_histograms/noise_histogram_g_right_down.npy', noise_histogram_g)
    np.save('noise_histograms/noise_histogram_bp_right_down.npy', noise_histogram_bp)
    np.save('noise_histograms/noise_histogram_rp_right_down.npy', noise_histogram_rp)
if screen_parts == 4 and histogram_parts == 100:
    area = ((max_ra - min_ra)/screen_parts * (max_dec - min_dec)/screen_parts)
    noise_histogram_g = np.array(noise_histogram_g) / area
    noise_histogram_bp = np.array(noise_histogram_bp) / area
    noise_histogram_rp = np.array(noise_histogram_rp) / area
    np.save('noise_histograms/noise_histogram_g_left_up.npy', noise_histogram_g)
    np.save('noise_histograms/noise_histogram_bp_left_up.npy', noise_histogram_bp)
    np.save('noise_histograms/noise_histogram_rp_left_up.npy', noise_histogram_rp)