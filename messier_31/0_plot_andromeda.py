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
# TODO ask file name instead of hardcoded

# filter results
results = results[results['ra_error'] > 0.1] # idk why this Gaia is so uncertain but smaller than 0.1 looks like noise
results = results[results['dec_error'] > 0.1]

# get useful data
ra = np.array(results['ra'])
dec = np.array(results['dec'])

phot_g_mean_flux = np.array(results['phot_g_mean_flux'])
phot_g_mean_mag = np.array(results['phot_g_mean_mag'])
# these are just the brightness, i need the start color "bp_rp" but many of the stars dont have those values :/
# i guess i have to work only with the brightness :(

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