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
results = results[~np.isnan(results['bp_rp'])]
results = results[~results['bp_rp'].mask]
results = results[~np.isnan(results['bp_g'])]
results = results[~results['bp_g'].mask]
results = results[~np.isnan(results['g_rp'])]
results = results[~results['g_rp'].mask]

# get useful data
ra = np.array(results['ra'])
dec = np.array(results['dec'])

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