import plotly.graph_objects as go
import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack

# import csv results
results = Table.read('andromeda_noiseless_rotated.csv', format='csv')

# prepare stars for histograme
threshhold = 0.1
screen_parts = 20
histogram_parts = 100
interesting_part = 1/4 # ratio of ra at which i want to compare the brightness of the stars with the opposite side

min_dec = np.min(results['dec'])
max_dec = np.max(results['dec'])
mid_dec = (max_dec + min_dec) / 2

results = results[(results['dec'] > mid_dec - threshhold) & (results['dec'] < mid_dec + threshhold)]

results.sort('ra')

# split X axis (ra) into screen_parts
min_ra = np.min(results['ra'])
max_ra = np.max(results['ra']) + 0.000001
bins_x = np.linspace(min_ra, max_ra, screen_parts+1)

split_results = []
for i in range(screen_parts):
    lower, upper = bins_x[i], bins_x[i+1]
    subset = results[(results['ra'] >= lower) & (results['ra'] < upper)]
    split_results.append(subset)

# create histogram inside each bin
histograms = []
min_brightness = np.min(results['phot_g_mean_mag'])
max_brightness = np.max(results['phot_g_mean_mag']) + 0.000001
for i in range(screen_parts):
    table = split_results[i]
    start_x = bins_x[i]
    end_x = bins_x[i+1]

    histogram = np.zeros(histogram_parts+1, dtype=int)
    for result in table:
        brightness = result['phot_g_mean_mag']
        brightness_ratio = (brightness - min_brightness) / (max_brightness - min_brightness)
        section = math.floor(brightness_ratio * histogram_parts)
        histogram[section] += 1
    sum_stars = np.sum(histogram)
    histogram = histogram / sum_stars
    histograms.append(histogram)

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



# plot the histograms
def num_to_color(num):
    if num < 0 or num > 1:
        return "gold"

    color_limit = 100  # from 0 to 255
    color_offset = 155  # color_offset + color_limit <= 255

    def to_hex(value):
        hex_value = hex(value)[2:]  # Convert to hex and remove "0x"
        return hex_value.zfill(2)  # Ensure two-digit format

    rgb = [
        math.floor(color_offset + color_limit / 2 * (math.sin((num      ) * 2 * math.pi) + 1)),
        math.floor(color_offset + color_limit / 2 * (math.sin((num + 1/3) * 2 * math.pi) + 1)),
        math.floor(color_offset + color_limit / 2 * (math.sin((num + 2/3) * 2 * math.pi) + 1))
    ]

    return "#" + "".join(to_hex(c) for c in rgb)

# Convert histograms into x, y, z data
x = []  # Histogram index
y = []  # Bin index
z = []  # Bin values (heights)

for hist_index, hist in enumerate(histograms):
    for bin_index, count in enumerate(hist):
        x.append(hist_index)  # Histogram number
        y.append(bin_index)   # Bin index in histogram
        z.append(count)       # Bin value (height)

# Create a 3D scatter plot using Plotly
fig = go.Figure(data=[go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        size=1,
        color=[num_to_color(i/2 / len(histograms)) for i in x],
        opacity=0.6
    ),
)])

# Labels and title
fig.update_layout(
    title="3D Plot of Star Positions",
    scene=dict(
        xaxis_title='Histogram Number (linear to ra)',
        yaxis_title='Brightness',
        zaxis_title='Number of stars normalized (num of stars / stars in histogram)',
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



# get interesting parts
histogram_num = math.floor(len(histograms) * interesting_part)
histogram1 = histograms[histogram_num]
histogram2 = histograms[len(histograms)-1 - histogram_num]
bins = np.arange(len(histogram1))

# Create a 2D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1,
    mode='markers',
    marker=dict(size=4, color='red'),
    line=dict(color='red', width=2),
    name="Histogram 1"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2,
    mode='markers',
    marker=dict(size=4, color='green'),
    line=dict(color='green', width=2),
    name="Histogram 2"
))

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
    ),
    yaxis=dict(
        showgrid=False,
        zeroline=False,
    ),
)

# Show the plot
fig.show()

# yep, using brightness is not eanough to tell difference between right and left size, means i cannot get any data from the speed of right or left sides.
# i need the color of the star, not just he brightness. but the database doesnt has colors for many stars.
# i wait until they improve data to include colors, or filter out only the stars that do have color data.

# TODO filter out the stars that only have their color (too lazy for now so i will leave this project until inspiration)



# bp_rp
# phot_g_mean_mag


# ok so i have colors data, i can do average and there will be difference, but the graphs look so random that even with that differents it looks not relaiable
# i think its bcz i randomly removed some of the noise, so now im also randomly left with some noise
# conclusion: i cannot randomly remove stuff, instead i should make histogram of the noise, histogram of the galaxy, and subtract. keep the big data big, no random shinanigans.

# also lets properly collect noise from both sides
# i will start over

# phot_g_n_obs
# phot_g_mean_flux
# phot_g_mean_flux_error
# phot_g_mean_flux_over_error
# phot_g_mean_mag
# phot_bp_mean_flux
# phot_bp_mean_flux_error
# phot_bp_mean_flux_over_error
# phot_bp_mean_mag
# phot_rp_mean_flux
# phot_rp_mean_flux_error
# phot_rp_mean_flux_over_error
# phot_rp_mean_mag
# phot_bp_rp_excess_factor
# bp_rp
# bp_g
# g_rp


# results = results[~np.isnan(results['bp_rp'])]
# results = results[~results['bp_rp'].mask]