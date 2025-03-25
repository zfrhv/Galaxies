import plotly.graph_objects as go
import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack

# import csv results
results = Table.read('andromeda_square_noiseless.csv', format='csv')

# load noise data
min_brightness = np.min(results['phot_g_mean_mag'])
max_brightness = np.max(results['phot_g_mean_mag']) + 0.000001
noise_histogram = np.load('noise_histogram_2.npy')

# filter out noise
noise_points_num = np.sum(noise_histogram)
random_x_s = np.random.uniform(low=np.min(results['ra']), high=np.max(results['ra']), size=noise_points_num)
random_y_s = np.random.uniform(low=np.min(results['dec']), high=np.max(results['dec']), size=noise_points_num)

results_array = results.as_array()
rows_to_delete = set()

j = 0
for i in range(len(noise_histogram)):
    points = noise_histogram[i]
    min_mag = min_brightness + (max_brightness - min_brightness) * (i / len(noise_histogram))
    max_mag = min_brightness + (max_brightness - min_brightness) * ((i + 1) / len(noise_histogram))
    mask = (results_array['phot_g_mean_mag'] >= min_mag) & (results_array['phot_g_mean_mag'] < max_mag)
    good_brightness = results_array[mask]

    # if edge of histogram with single point then multiplied by 9 then just remove those points
    if len(good_brightness) < points:
        points = len(good_brightness)

    for point in range(points):
        random_x = random_x_s[j]
        random_y = random_y_s[j]

        distances = (good_brightness['ra'] - random_x) ** 2 + (good_brightness['dec'] - random_y) ** 2

        # Find closest star thats not already deleted
        indexes = np.argsort(distances)
        real_index = None
        found = False
        for index in indexes:
            real_index = np.where(mask)[0][index]
            if real_index not in rows_to_delete:
                found = True
                break

        if found:
            rows_to_delete.add(real_index)
        else:
            print('error: couldnt find fitting star')

        j += 1
        if j % 100 == 0:
            print(j/noise_points_num*100, '%')

    # clean trash if big
    if len(rows_to_delete) >= 100:
        print('cleaning the trash')
        mask_2 = np.ones(len(results), dtype=bool)
        mask_2[list(rows_to_delete)] = False
        results = results[mask_2]
        results_array = results.as_array()
        rows_to_delete = set()

# Step 3: Efficiently remove all selected rows in one go
mask_2 = np.ones(len(results), dtype=bool)
mask_2[list(rows_to_delete)] = False
results = results[mask_2]

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

# save the table
results.write('andromeda_square_noiseless_2.csv', format='csv', overwrite=True)