import plotly.graph_objects as go
import numpy as np
import math
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import gaussian_filter

# import csv results
results = Table.read('andromeda_rotated.csv', format='csv')

# prepare stars for histograme
threshhold = 0.1
screen_parts = 20
histogram_parts = 1000
interesting_part = 5 # which part from left to compare with with its mirror part on right indexes are 0-19

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

    histogram_bp_rp = np.zeros(histogram_parts+1, dtype=float)
    histogram_bp_g = np.zeros(histogram_parts+1, dtype=float)
    histogram_g_rp = np.zeros(histogram_parts+1, dtype=float)
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


# load noise data
noise_histogram_bp_rp_left_up = np.load('noise_histograms/noise_histogram_bp_rp_left_up.npy')
noise_histogram_bp_g_left_up = np.load('noise_histograms/noise_histogram_bp_g_left_up.npy')
noise_histogram_g_rp_left_up = np.load('noise_histograms/noise_histogram_g_rp_left_up.npy')
noise_area_left_up = np.load('noise_histograms/noise_area_left_up.npy')


# seperate histograms into left, avg, right
hist_parts = math.floor(screen_parts/2)
noise_area = noise_area_left_up
bp_rp = {
    "left_hists": histograms_bp_rp[:hist_parts],
    "right_hists": histograms_bp_rp[hist_parts:],
    "noise_hist": noise_histogram_bp_rp_left_up
}
bp_rp["avg_hists"] = [(h1 + h2) / 2 for h1, h2 in zip(bp_rp["left_hists"], bp_rp["right_hists"])]

bp_g = {
    "left_hists": histograms_bp_g[:hist_parts],
    "right_hists": histograms_bp_g[hist_parts:],
    "noise_hist": noise_histogram_bp_g_left_up
}
bp_g["avg_hists"] = [(h1 + h2) / 2 for h1, h2 in zip(bp_g["left_hists"], bp_g["right_hists"])]

g_rp = {
    "left_hists": histograms_g_rp[:hist_parts],
    "right_hists": histograms_g_rp[hist_parts:],
    "noise_hist": noise_histogram_g_rp_left_up
}
g_rp["avg_hists"] = [(h1 + h2) / 2 for h1, h2 in zip(g_rp["left_hists"], g_rp["right_hists"])]


# smooth histograms
def smooth_histogram(hist):
    hist[:] = gaussian_filter1d(hist, sigma=10)


for color_shade in [bp_rp, bp_g, g_rp]:
    smooth_histogram(color_shade["noise_hist"])

# smooth between the histograms
def smooth_between(hists):
    hists[:] = gaussian_filter(hists, sigma=(1.5, 10))
    # hists[:] = gaussian_filter(hists, sigma=(0, 10))

for color_shade in [bp_rp, bp_g, g_rp]:
    smooth_between(color_shade["left_hists"])
    smooth_between(color_shade["right_hists"])
    smooth_between(color_shade["avg_hists"])

# remove the noise
new_area = (max_ra - min_ra)/screen_parts * (threshhold*2)
area_ratio = new_area / noise_area

# TODO get this back?
# for color_shade in [bp_rp, bp_g, g_rp]:
#     color_shade["left_hists"]  = color_shade["left_hists"]  - color_shade["noise_hist"] * area_ratio
#     color_shade["right_hists"] = color_shade["right_hists"] - color_shade["noise_hist"] * area_ratio
#     color_shade["avg_hists"]   = color_shade["avg_hists"]   - color_shade["noise_hist"] * area_ratio

# # straight the amount of stars, to get avg colors per 1 star
# for color_shade in [bp_rp, bp_g, g_rp]:
#     for hists in [color_shade["left_hists"], color_shade["right_hists"], color_shade["avg_hists"]]:
#         for i in range(len(hists)):
#             hists[i] = hists[i] / np.sum(hists[i])


















# get useful data
ra = np.array(original_results['ra'])
dec = np.array(original_results['dec'])
colors = np.array(original_results['bp_rp'])

# Create a 2D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scattergl(
    x=ra, y=dec,
    mode='markers',
    marker=dict(
        size=1,
        color=colors,
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
fig = go.Figure()

for color_shade, color in zip([bp_rp, bp_g, g_rp], ["green", "lightblue", "red"]):
    # Convert histograms into x, y, z data
    x = [] # Histogram index
    y = [] # Bin index
    z = [] # Bin values (heights)
    for hist_index, hist in enumerate(color_shade["left_hists"]):
        for bin_index, count in enumerate(hist):
            x.append(hist_index)  # Histogram number
            y.append(bin_index)   # Bin index in histogram
            z.append(count)       # Bin value (height)

    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=1,
            color=color,
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
fig = go.Figure()

for color_shade, name, color in zip([bp_rp, bp_g, g_rp], ["bp_rp", "bp_g", "g_rp"], ["green", "lightblue", "red"]):
    left_hist = color_shade["left_hists"][interesting_part]
    right_hist = color_shade["right_hists"][interesting_part]
    avg_hist = color_shade["avg_hists"][interesting_part]
    # TODO check this minus
    # left_hist = left_hist - avg_hist
    # right_hist = right_hist - avg_hist
    bins = np.arange(len(left_hist))

    print('speed',name,':', sum(val * i for i, val in enumerate(right_hist)) / 6.499625227109498)

    fig.add_trace(go.Scattergl(
        x=bins, y=left_hist,
        mode='lines+markers',
        marker=dict(size=1, color=color),
        name=name+" left"
    ))

    fig.add_trace(go.Scattergl(
        x=bins, y=right_hist,
        mode='lines+markers',
        marker=dict(size=1, color=color),
        name=name+" right"
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