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

# rotate the results
degree = math.radians(-45)

min_ra = np.min(results['ra'])
max_ra = np.max(results['ra'])
min_dec = np.min(results['dec'])
max_dec = np.max(results['dec'])
mid_ra = (max_ra + min_ra)/2
mid_dec = (max_dec + min_dec)/2

for result in results:
    x = result['ra']
    y = result['dec']
    cos_theta = math.cos(degree)
    sin_theta = math.sin(degree)
    x_new = x * cos_theta - y * sin_theta
    y_new = x * sin_theta + y * cos_theta
    result['ra'] = x_new
    result['dec'] = y_new

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
results.write('andromeda_rotated.csv', format='csv', overwrite=True)