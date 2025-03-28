import plotly.graph_objects as go
import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack

# import csv results
results = Table.read('andromeda_rotated.csv', format='csv')

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
histograms_bp_rp = []
histograms_bp_g = []
histograms_g_rp = []
min_bp_rp = np.min(results['bp_rp'])
max_bp_rp = np.max(results['bp_rp']) + 0.000001
min_bp_g = np.min(results['bp_g'])
max_bp_g = np.max(results['bp_g']) + 0.000001
min_g_rp = np.min(results['g_rp'])
max_g_rp = np.max(results['g_rp']) + 0.000001
for i in range(screen_parts):
    table = split_results[i]
    start_x = bins_x[i]
    end_x = bins_x[i+1]

    histogram_bp_rp = np.zeros(histogram_parts+1, dtype=int)
    histogram_bp_g = np.zeros(histogram_parts+1, dtype=int)
    histogram_g_rp = np.zeros(histogram_parts+1, dtype=int)
    for result in table:
        bp_rp = result['bp_rp']
        bp_rp_ratio = (bp_rp - min_bp_rp) / (max_bp_rp - min_bp_rp)
        section_bp_rp = math.floor(bp_rp_ratio * histogram_parts)
        histogram_bp_rp[section_bp_rp] += 1

        bp_g = result['bp_g']
        bp_g_ratio = (bp_g - min_bp_g) / (max_bp_g - min_bp_g)
        section_bp_g = math.floor(bp_g_ratio * histogram_parts)
        histogram_bp_g[section_bp_g] += 1

        g_rp = result['g_rp']
        g_rp_ratio = (g_rp - min_g_rp) / (max_g_rp - min_g_rp)
        section_g_rp = math.floor(g_rp_ratio * histogram_parts)
        histogram_g_rp[section_g_rp] += 1

    sum_stars_bp_rp = np.sum(histogram_bp_rp)
    histogram_bp_rp = histogram_bp_rp / sum_stars_bp_rp
    histograms_bp_rp.append(histogram_bp_rp)

    sum_stars_bp_g = np.sum(histogram_bp_g)
    histogram_bp_g = histogram_bp_g / sum_stars_bp_g
    histograms_bp_g.append(histogram_bp_g)

    sum_stars_g_rp = np.sum(histogram_g_rp)
    histogram_g_rp = histogram_g_rp / sum_stars_g_rp
    histograms_g_rp.append(histogram_g_rp)


# remove the noise
noise_histogram_bp_rp_right_down = np.load('noise_histograms/noise_histogram_bp_rp_right_down.npy')
noise_histogram_bp_g_right_down = np.load('noise_histograms/noise_histogram_bp_g_right_down.npy')
noise_histogram_g_rp_right_down = np.load('noise_histograms/noise_histogram_g_rp_right_down.npy')
noise_histogram_bp_rp_left_up = np.load('noise_histograms/noise_histogram_bp_rp_left_up.npy')
noise_histogram_bp_g_left_up = np.load('noise_histograms/noise_histogram_bp_g_left_up.npy')
noise_histogram_g_rp_left_up = np.load('noise_histograms/noise_histogram_g_rp_left_up.npy')

# noise per area
area = (max_ra - min_ra)/screen_parts * (threshhold*2)
noise_histogram_bp_rp = (noise_histogram_bp_rp_right_down + noise_histogram_bp_rp_left_up) / 2 * area
noise_histogram_bp_g = (noise_histogram_bp_g_right_down + noise_histogram_bp_g_left_up) / 2 * area
noise_histogram_g_rp = (noise_histogram_g_rp_right_down + noise_histogram_g_rp_left_up) / 2 * area

for i in range(screen_parts):
    histograms_bp_rp[i] = histograms_bp_rp[i] - noise_histogram_bp_rp
    histograms_bp_g[i] = histograms_bp_g[i] - noise_histogram_bp_g
    histograms_g_rp[i] = histograms_g_rp[i] - noise_histogram_g_rp
















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






















# plot the histograms
# Convert histograms into x, y, z data
bp_rp_x = []  # Histogram index
bp_rp_y = []  # Bin index
bp_rp_z = []  # Bin values (heights)
bp_g_x = []
bp_g_y = []
bp_g_z = []
g_rp_x = []
g_rp_y = []
g_rp_z = []

for hist_index, hist in enumerate(histograms_bp_rp):
    for bin_index, count in enumerate(hist):
        bp_rp_x.append(hist_index)  # Histogram number
        bp_rp_y.append(bin_index)   # Bin index in histogram
        bp_rp_z.append(count)       # Bin value (height)
for hist_index, hist in enumerate(histograms_bp_g):
    for bin_index, count in enumerate(hist):
        bp_g_x.append(hist_index)
        bp_g_y.append(bin_index)
        bp_g_z.append(count)
for hist_index, hist in enumerate(histograms_g_rp):
    for bin_index, count in enumerate(hist):
        g_rp_x.append(hist_index)
        g_rp_y.append(bin_index)
        g_rp_z.append(count)

# Create a 3D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=bp_rp_x,
    y=bp_rp_y,
    z=bp_rp_z,
    mode='markers',
    marker=dict(
        size=1,
        color='green',
        opacity=0.6
    ),
))
fig.add_trace(go.Scatter3d(
    x=bp_g_x,
    y=bp_g_y,
    z=bp_g_z,
    mode='markers',
    marker=dict(
        size=1,
        color='lightblue',
        opacity=0.6
    ),
))
fig.add_trace(go.Scatter3d(
    x=g_rp_x,
    y=g_rp_y,
    z=g_rp_z,
    mode='markers',
    marker=dict(
        size=1,
        color='red',
        opacity=0.6
    ),
))

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

# Create a 2D scatter plot using Plotly
fig = go.Figure()

histogram_num = math.floor(len(histograms_bp_rp) * interesting_part)
histogram1_bp_rp = histograms_bp_rp[histogram_num]
histogram2_bp_rp = histograms_bp_rp[len(histograms_bp_rp)-1 - histogram_num]
bins = np.arange(len(histogram1_bp_rp))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_bp_rp,
    mode='lines+markers',
    marker=dict(size=4, color='red'),
    line=dict(color='green', width=2),
    name="green left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_bp_rp,
    mode='lines+markers',
    marker=dict(size=4, color='lightblue'),
    line=dict(color='green', width=2),
    name="green right"
))

histogram_num = math.floor(len(histograms_bp_g) * interesting_part)
histogram1_bp_g = histograms_bp_g[histogram_num]
histogram2_bp_g = histograms_bp_g[len(histograms_bp_g)-1 - histogram_num]
bins = np.arange(len(histogram1_bp_g))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_bp_g,
    mode='lines+markers',
    marker=dict(size=4, color='red'),
    line=dict(color='lightblue', width=2),
    name="blue left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_bp_g,
    mode='lines+markers',
    marker=dict(size=4, color='green'),
    line=dict(color='lightblue', width=2),
    name="blue right"
))

histogram_num = math.floor(len(histograms_g_rp) * interesting_part)
histogram1_g_rp = histograms_g_rp[histogram_num]
histogram2_g_rp = histograms_g_rp[len(histograms_g_rp)-1 - histogram_num]
bins = np.arange(len(histogram1_g_rp))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_g_rp,
    mode='lines+markers',
    marker=dict(size=4, color='green'),
    line=dict(color='red', width=2),
    name="red left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_g_rp,
    mode='lines+markers',
    marker=dict(size=4, color='lightblue'),
    line=dict(color='red', width=2),
    name="red right"
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



# ah nvm they are also too close
# and the left vs right parts also keep switching values so i dont think its relayable
# i dont think i can do anything with this data, there is no color difference between left and right sides
# so i cant tells any speed differences between left and right
