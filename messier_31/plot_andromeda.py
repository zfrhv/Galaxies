import plotly.graph_objects as go
import numpy as np
from astropy.io import fits
from astropy.table import Table

# import results
hdul = fits.open('andromeda_data.fits')
data = hdul[1].data
results = Table(data)
hdul.close()

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

distance = 2.537 * 10**6 # light years

# Convert spherical coordinates (RA, DEC, Distance) to Cartesian coordinates (x, y, z)
ra_rad = np.radians(ra)
dec_rad = np.radians(dec)
x = distance * np.cos(dec_rad) * np.cos(ra_rad)
y = distance * np.cos(dec_rad) * np.sin(ra_rad)
z = distance * np.sin(dec_rad)

# Create a 3D scatter plot using Plotly
fig = go.Figure(data=[go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=1,
        color='white',
        colorscale='YlGnBu',
        opacity=0.6
    ),
)])

# Labels and title
fig.update_layout(
    title="3D Plot of Star Positions",
    scene=dict(
        xaxis_title='X [ly]',
        yaxis_title='Y [ly]',
        zaxis_title='Z [ly]',
        bgcolor='black',
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showbackground=False,
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showbackground=False,
        ),
        zaxis=dict(
            showgrid=False,
            zeroline=False,
            showbackground=False,
        ),
    ),
    plot_bgcolor='black',
    paper_bgcolor='black',
    font=dict(color='white'),
)

# Show the plot
fig.show()