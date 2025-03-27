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
results = results[~np.isnan(results['phot_g_mean_flux'])]
results = results[~np.isnan(results['phot_g_mean_mag'])]
results = results[~results['phot_g_mean_flux'].mask]
results = results[~results['phot_g_mean_mag'].mask]

# sort stars for histograme
results.sort('ra')

screen_parts = 3
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
noise_histogram = None
min_brightness = np.min(results['phot_g_mean_mag'])
max_brightness = np.max(results['phot_g_mean_mag']) + 0.000001
for i in range(screen_parts):
    for j in range(screen_parts):
        table = split_results[i][j]
        start_x = bins_x[i]
        end_x = bins_x[i+1]
        start_y = bins_y[j]
        end_y = bins_y[j+1]

        histogram = np.zeros(histogram_parts+1, dtype=int)
        section_size = (max_ra - min_ra) / screen_parts / histogram_parts
        for result in table:
            brightness = result['phot_g_mean_mag']
            brightness_ratio = (brightness - min_brightness) / (max_brightness - min_brightness)
            section = math.floor(brightness_ratio * histogram_parts)
            offset_x = section * section_size
            offset_y = histogram[section] * (section_size * histogram_height)
            histogram[section] += 1

            result['ra'] = start_x + offset_x
            result['dec'] = start_y + offset_y
        # use the right bottom histogram for noise inspection
        if i == screen_parts-1 and j == 0:
            noise_histogram = histogram



# merge the results back into one table
flattened_results = []
for sublist in split_results:
    for table in sublist:
        flattened_results.append(table)
split_results = flattened_results

results = vstack(split_results)

# get useful data
ra = np.array(results['ra'])
dec = np.array(results['dec'])
phot_g_mean_flux = np.array(results['phot_g_mean_flux'])
phot_g_mean_mag = np.array(results['phot_g_mean_mag'])

# Create a 2D scatter plot using Plotly
fig = go.Figure(data=[go.Scattergl(
    x=ra,
    y=dec,
    mode='markers',
    marker=dict(
        size=1,
        color=phot_g_mean_mag,
        colorscale='hot',
        opacity=1
    ),
)])

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
if screen_parts == 3:
    noise_histogram = np.array(noise_histogram) * 9 # multiply by 3*3 cuz this histogram is only 1/9 of the screen
    np.save('noise_histogram.npy', noise_histogram)

# TODO maybe use phot_g_mean_flux instead?