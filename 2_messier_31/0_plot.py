import plotly.graph_objects as go
import numpy as np
from astropy.io import fits
from astropy.table import Table
import tkinter as tk
from tkinter import filedialog

# import csv results
filename = filedialog.askopenfilename(filetypes=[("CSV files", ".csv")])
results = Table.read(filename, format='csv')

# filter results
results = results[results['ra_error'] > 0.1] # idk why this Gaia is so uncertain but smaller than 0.1 looks like noise
results = results[results['dec_error'] > 0.1]
results = results[~np.isnan(results['phot_g_mean_flux'])]
results = results[~results['phot_g_mean_flux'].mask]
results = results[~np.isnan(results['phot_bp_mean_flux'])]
results = results[~results['phot_bp_mean_flux'].mask]
results = results[~np.isnan(results['phot_rp_mean_flux'])]
results = results[~results['phot_rp_mean_flux'].mask]
results = results[~np.isnan(results['phot_g_mean_flux_error'])]
results = results[~results['phot_g_mean_flux_error'].mask]
results = results[~np.isnan(results['phot_bp_mean_flux_error'])]
results = results[~results['phot_bp_mean_flux_error'].mask]
results = results[~np.isnan(results['phot_rp_mean_flux_error'])]
results = results[~results['phot_rp_mean_flux_error'].mask]

# get useful data
ra = np.array(results['ra'])
dec = np.array(results['dec'])
phot_g_mean_flux = np.array(results['phot_g_mean_flux'])
phot_bp_mean_flux = np.array(results['phot_bp_mean_flux'])
phot_rp_mean_flux = np.array(results['phot_rp_mean_flux'])
phot_g_mean_flux_error = np.array(results['phot_g_mean_flux_error'])
phot_bp_mean_flux_error = np.array(results['phot_bp_mean_flux_error'])
phot_rp_mean_flux_error = np.array(results['phot_rp_mean_flux_error'])

# Create a 2D scatter plot using Plotly
fig = go.Figure(data=[go.Scattergl(
    x=ra,
    y=dec,
    mode='markers',
    marker=dict(
        size=1,
        color='white',
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