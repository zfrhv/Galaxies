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
histogram_parts = 1000
interesting_part = 3 # which part from left to compare with with its mirror part on right indexes are 0-19

min_dec = np.min(results['dec'])
max_dec = np.max(results['dec'])
mid_dec = (max_dec + min_dec) / 2

original_results = results
results = results[(results['dec'] > mid_dec - threshhold) & (results['dec'] < mid_dec + threshhold)]

results.sort('ra')

# split X axis (ra) into screen_parts
padding = threshhold + 0.01 # bcz picture rotated by 45.5 degress
min_ra = np.min(results['ra']) + padding
max_ra = np.max(results['ra']) - padding + 0.000001
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
        section = math.floor(bp_rp_ratio * histogram_parts)
        histogram_bp_rp[section] += 1

        bp_g = result['bp_g']
        bp_g_ratio = (bp_g - min_bp_g) / (max_bp_g - min_bp_g)
        section = math.floor(bp_g_ratio * histogram_parts)
        histogram_bp_g[section] += 1

        g_rp = result['g_rp']
        g_rp_ratio = (g_rp - min_g_rp) / (max_g_rp - min_g_rp)
        section = math.floor(g_rp_ratio * histogram_parts)
        histogram_g_rp[section] += 1

    histograms_bp_rp.append(histogram_bp_rp)
    histograms_bp_g.append(histogram_bp_g)
    histograms_g_rp.append(histogram_g_rp)


# remove the noise
noise_histogram_bp_rp_left_up = np.load('noise_histograms/noise_histogram_bp_rp_left_up.npy')
noise_histogram_bp_g_left_up = np.load('noise_histograms/noise_histogram_bp_g_left_up.npy')
noise_histogram_g_rp_left_up = np.load('noise_histograms/noise_histogram_g_rp_left_up.npy')

# noise per area
area = (max_ra - min_ra)/screen_parts * (threshhold*2)
noise_histogram_bp_rp = noise_histogram_bp_rp_left_up * area
noise_histogram_bp_g = noise_histogram_bp_g_left_up * area
noise_histogram_g_rp = noise_histogram_g_rp_left_up * area

for i in range(screen_parts):
    histograms_bp_rp[i] = histograms_bp_rp[i] - noise_histogram_bp_rp
    histograms_bp_g[i] = histograms_bp_g[i] - noise_histogram_bp_g
    histograms_g_rp[i] = histograms_g_rp[i] - noise_histogram_g_rp

# straight the amount of stars, to get avg colors per 1 star
for i in range(screen_parts):
    histograms_bp_rp[i] = histograms_bp_rp[i] / np.sum(histograms_bp_rp[i])
    histograms_bp_g[i] = histograms_bp_g[i] / np.sum(histograms_bp_g[i])
    histograms_g_rp[i] = histograms_g_rp[i] / np.sum(histograms_g_rp[i])


# seperate histograms into left, avg, right
hist_parts = math.floor(screen_parts/2)
left_hist_bp_rp = histograms_bp_rp[:hist_parts]
right_hist_bp_rp = histograms_bp_rp[hist_parts:]
avg_hist_bp_rp = [(h1 + h2) / 2 for h1, h2 in zip(left_hist_bp_rp, right_hist_bp_rp)]

left_hist_bp_g = histograms_bp_g[:hist_parts]
right_hist_bp_g = histograms_bp_g[hist_parts:]
avg_hist_bp_g = [(h1 + h2) / 2 for h1, h2 in zip(left_hist_bp_g, right_hist_bp_g)]

left_hist_g_rp = histograms_g_rp[:hist_parts]
right_hist_g_rp = histograms_g_rp[hist_parts:]
avg_hist_g_rp = [(h1 + h2) / 2 for h1, h2 in zip(left_hist_g_rp, right_hist_g_rp)]

# smooth each histogram
def smooth_histogram(hist):
    # take point befor and point after, not avg but also check the speeds in range..
    # dont forget that when i modify value the next itteration will use it as anchor
    print()

for hist in left_hist_bp_rp:
    smooth_histogram(hist)
for hist in right_hist_bp_rp:
    smooth_histogram(hist)
for hist in avg_hist_bp_rp:
    smooth_histogram(hist)

for hist in left_hist_bp_g:
    smooth_histogram(hist)
for hist in right_hist_bp_g:
    smooth_histogram(hist)
for hist in avg_hist_bp_g:
    smooth_histogram(hist)

for hist in left_hist_g_rp:
    smooth_histogram(hist)
for hist in right_hist_g_rp:
    smooth_histogram(hist)
for hist in avg_hist_g_rp:
    smooth_histogram(hist)

# smooth between the histograms
def smooth_between(hists):
    print()

smooth_between(left_hist_bp_rp)
smooth_between(right_hist_bp_rp)
smooth_between(avg_hist_bp_rp)

smooth_between(left_hist_bp_g)
smooth_between(right_hist_bp_g)
smooth_between(avg_hist_bp_g)

smooth_between(left_hist_g_rp)
smooth_between(right_hist_g_rp)
smooth_between(avg_hist_g_rp)

















# get useful data
ra = np.array(original_results['ra'])
dec = np.array(original_results['dec'])
bp_rp = np.array(original_results['bp_rp'])

# Create a 2D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scattergl(
    x=ra,
    y=dec,
    mode='markers',
    marker=dict(
        size=1,
        color=bp_rp,
        colorscale='hot',
        opacity=1
    ),
))

# Add borders
borders_x = [min_ra, max_ra, max_ra, min_ra, min_ra]
borders_y = [mid_dec + threshhold, mid_dec + threshhold, mid_dec - threshhold, mid_dec - threshhold, mid_dec + threshhold]
for i in range(1, screen_parts):
    new_x = (max_ra - min_ra)/screen_parts * i + min_ra
    borders_x.extend([new_x, new_x, new_x])
    borders_y.extend([mid_dec + threshhold, mid_dec - threshhold, mid_dec + threshhold])
fig.add_trace(go.Scattergl(
    x=borders_x,
    y=borders_y,
    mode='lines',
    line=dict(color='white', width=1),
    name="Square"
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

for hist_index, hist in enumerate(left_hist_bp_rp):
    for bin_index, count in enumerate(hist):
        bp_rp_x.append(hist_index)  # Histogram number
        bp_rp_y.append(bin_index)   # Bin index in histogram
        bp_rp_z.append(count)       # Bin value (height)
for hist_index, hist in enumerate(left_hist_bp_g):
    for bin_index, count in enumerate(hist):
        bp_g_x.append(hist_index)
        bp_g_y.append(bin_index)
        bp_g_z.append(count)
for hist_index, hist in enumerate(left_hist_g_rp):
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

left_hist = left_hist_bp_rp[interesting_part]
right_hist = right_hist_bp_rp[interesting_part]
avg_hist = avg_hist_bp_rp[interesting_part]
left_hist = left_hist - avg_hist
right_hist = right_hist - avg_hist
bins = np.arange(len(left_hist))

print('speed by bp_rp: ', sum(val * i for i, val in enumerate(right_hist)) / 6.499625227109498)

fig.add_trace(go.Scattergl(
    x=bins,
    y=left_hist,
    mode='lines+markers',
    marker=dict(size=1, color='green'),
    name="bp_rp left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=right_hist,
    mode='lines+markers',
    marker=dict(size=1, color='lightgreen'),
    name="bp_rp right"
))

left_hist = left_hist_bp_g[interesting_part]
right_hist = right_hist_bp_g[interesting_part]
avg_hist = avg_hist_bp_g[interesting_part]
left_hist = left_hist - avg_hist
right_hist = right_hist - avg_hist
bins = np.arange(len(left_hist))

print('speed by bp_g: ', sum(val * i for i, val in enumerate(right_hist)) / 2.3683722618236382)

fig.add_trace(go.Scattergl(
    x=bins,
    y=left_hist,
    mode='lines+markers',
    marker=dict(size=1, color='blue'),
    name="bp_g left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=right_hist,
    mode='lines+markers',
    marker=dict(size=1, color='lightblue'),
    name="bp_g right"
))

left_hist = left_hist_g_rp[interesting_part]
right_hist = right_hist_g_rp[interesting_part]
avg_hist = avg_hist_g_rp[interesting_part]
left_hist = left_hist - avg_hist
right_hist = right_hist - avg_hist
bins = np.arange(len(left_hist))

print('speed by g_rp: ', sum(val * i for i, val in enumerate(right_hist)) / 4.13610394421988)

fig.add_trace(go.Scattergl(
    x=bins,
    y=left_hist,
    mode='lines+markers',
    marker=dict(size=1, color='red'),
    name="g_rp left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=right_hist,
    mode='lines+markers',
    marker=dict(size=1, color='pink'),
    name="g_rp right"
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



# TODO
# make smooth_histogram() and smooth all histograms
# smooth between the left histograms same for avg and right
# do left-avg vs right-avg and see the speeds