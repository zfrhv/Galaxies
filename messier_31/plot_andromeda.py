import plotly.graph_objects as go
import numpy as np
from astropy.io import fits
from astropy.table import Table

hdul = fits.open('andromeda_data.fits')
data = hdul[1].data
results = Table(data)
hdul.close()

ra = np.array(results['ra'])
dec = np.array(results['dec'])
parallax = np.array(results['parallax'])
flux = np.array(results['phot_g_mean_flux'])

min_flux = np.min(flux)
max_flux = np.max(flux)
brightness = (flux - (min_flux / 2)) / (max_flux - min_flux)

# Convert RA and DEC to radians
ra_rad = np.radians(ra)
dec_rad = np.radians(dec)

distance = 2.537 # million light years

# Convert spherical coordinates (RA, DEC, Distance) to Cartesian coordinates (x, y, z)
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
        color=brightness,
        colorscale='YlGnBu',
        opacity=1
    ),
)])

# Labels and title
fig.update_layout(
    title="3D Plot of Star Positions",
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
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