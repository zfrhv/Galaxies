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

# sort stars for histograme
results.sort('ra')
bin_size = 0.2
histogram_parts = 20

histogram_height = histogram_parts * 1/bin_size * 1/10000

# split X axis (ra) into small parts size of bin_size
min_ra = np.min(results['ra'])
max_ra = np.max(results['ra'])
bins_x = np.arange(min_ra, max_ra+bin_size, bin_size)

split_results = []
for i in range(len(bins_x) - 1):
    lower, upper = bins_x[i], bins_x[i+1]
    subset = results[(results['ra'] >= lower) & (results['ra'] < upper)]
    split_results.append(subset)

# split Y axis (dec) into small parts size of bin_size
min_dec = np.min(results['dec'])
max_dec = np.max(results['dec'])
bins_y = np.arange(min_dec, max_dec+bin_size, bin_size)

for i in range(len(split_results)):
    sub_results = split_results[i]
    split_sub_results = []
    for j in range(len(bins_y) - 1):
        lower, upper = bins_y[j], bins_y[j+1]
        subset = sub_results[(sub_results['dec'] >= lower) & (sub_results['dec'] < upper)]
        split_sub_results.append(subset)
    split_results[i] = split_sub_results

# create histogram inside each bin
min_brightness = np.min(results['phot_g_mean_mag'])
max_brightness = np.max(results['phot_g_mean_mag'])
for i in range(len(bins_x) - 1):
    for j in range(len(bins_y) - 1):
        table = split_results[i][j]
        start_x = bins_x[i]
        end_x = bins_x[i+1]
        start_y = bins_y[j]
        end_y = bins_y[j+1]

        histogram = np.zeros(histogram_parts+1)
        section_size = bin_size / histogram_parts
        for result in table:
            brightness = result['phot_g_mean_mag']
            if np.ma.is_masked(brightness):
                print('Warning: some brightness is masked')
                continue
            brightness_ratio = (brightness - min_brightness) / (max_brightness - min_brightness)
            section = math.floor(brightness_ratio * histogram_parts)
            offset_x = section * section_size
            offset_y = histogram[section] * (section_size * histogram_height)
            histogram[section] += 1

            result['ra'] = start_x + offset_x
            result['dec'] = start_y + offset_y



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