import plotly.graph_objects as go
import numpy as np
from astropy.io import fits
from astropy.table import Table

# # import fits results
# hdul = fits.open('andromeda_cone_data.fits')
# data = hdul[1].data
# results = Table(data)
# hdul.close()

# import csv results
results = Table.read('andromeda_square_data.csv', format='csv')

# filter results
results = results[results['ra_error'] > 0.1] # idk why this Gaia is so uncertain but smaller than 0.1 looks like noise
results = results[results['dec_error'] > 0.1]

# get useful data
ra = np.array(results['ra'])
dec = np.array(results['dec'])

phot_g_mean_flux = np.array(results['phot_g_mean_flux'])
phot_g_mean_mag = np.array(results['phot_g_mean_mag'])

# astrometric_gof_al : 14.049725
# astrometric_chi2_al : 348.53717
# astrometric_excess_noise : 25.398783
# astrometric_excess_noise_sig : 21.824903
# astrometric_sigma5d_max : 948.63873
# ipd_gof_harmonic_amplitude : 1.9720112
# ipd_gof_harmonic_phase : 50.664013

# # rotate to make it easier
# theta = np.radians(-37)
# ra = ra * np.cos(theta) - dec * np.sin(theta)
# dec = ra * np.sin(theta) + dec * np.cos(theta)

# Create a 2D scatter plot using Plotly
fig = go.Figure(data=[go.Scatter(
    x=ra,
    y=dec,
    mode='markers',
    marker=dict(
        size=1,
        color='white',
        opacity=0.6
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