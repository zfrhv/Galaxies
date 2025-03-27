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
histograms_g = []
histograms_bp = []
histograms_rp = []
min_phot_g = np.min(results['phot_g_mean_mag'])
max_phot_g = np.max(results['phot_g_mean_mag']) + 0.000001
min_phot_bp = np.min(results['phot_bp_mean_mag'])
max_phot_bp = np.max(results['phot_bp_mean_mag']) + 0.000001
min_phot_rp = np.min(results['phot_rp_mean_mag'])
max_phot_rp = np.max(results['phot_rp_mean_mag']) + 0.000001
for i in range(screen_parts):
    table = split_results[i]
    start_x = bins_x[i]
    end_x = bins_x[i+1]

    histogram_g = np.zeros(histogram_parts+1, dtype=int)
    histogram_bp = np.zeros(histogram_parts+1, dtype=int)
    histogram_rp = np.zeros(histogram_parts+1, dtype=int)
    for result in table:
        phot_g = result['phot_g_mean_mag']
        phot_g_ratio = (phot_g - min_phot_g) / (max_phot_g - min_phot_g)
        section_g = math.floor(phot_g_ratio * histogram_parts)
        histogram_g[section_g] += 1

        phot_bp = result['phot_bp_mean_mag']
        phot_bp_ratio = (phot_bp - min_phot_bp) / (max_phot_bp - min_phot_bp)
        section_bp = math.floor(phot_bp_ratio * histogram_parts)
        histogram_bp[section_bp] += 1

        phot_rp = result['phot_rp_mean_mag']
        phot_rp_ratio = (phot_rp - min_phot_rp) / (max_phot_rp - min_phot_rp)
        section_rp = math.floor(phot_rp_ratio * histogram_parts)
        histogram_rp[section_rp] += 1

    sum_stars_g = np.sum(histogram_g)
    histogram_g = histogram_g / sum_stars_g
    histograms_g.append(histogram_g)

    sum_stars_bp = np.sum(histogram_bp)
    histogram_bp = histogram_bp / sum_stars_bp
    histograms_bp.append(histogram_bp)

    sum_stars_rp = np.sum(histogram_rp)
    histogram_rp = histogram_rp / sum_stars_rp
    histograms_rp.append(histogram_rp)


# remove the noise
noise_histogram_g_right_down = np.load('noise_histograms/noise_histogram_g_right_down.npy')
noise_histogram_bp_right_down = np.load('noise_histograms/noise_histogram_bp_right_down.npy')
noise_histogram_rp_right_down = np.load('noise_histograms/noise_histogram_rp_right_down.npy')
noise_histogram_g_left_up = np.load('noise_histograms/noise_histogram_g_left_up.npy')
noise_histogram_bp_left_up = np.load('noise_histograms/noise_histogram_bp_left_up.npy')
noise_histogram_rp_left_up = np.load('noise_histograms/noise_histogram_rp_left_up.npy')

# noise per area
area = (max_ra - min_ra)/screen_parts * (threshhold*2)
noise_histogram_g = (noise_histogram_g_right_down + noise_histogram_g_left_up) / 2 * area
noise_histogram_bp = (noise_histogram_bp_right_down + noise_histogram_bp_left_up) / 2 * area
noise_histogram_rp = (noise_histogram_rp_right_down + noise_histogram_rp_left_up) / 2 * area

for i in range(screen_parts):
    histograms_g[i] = histograms_g[i] - noise_histogram_g
    histograms_bp[i] = histograms_bp[i] - noise_histogram_bp
    histograms_rp[i] = histograms_rp[i] - noise_histogram_rp
















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
phot_g_x = []  # Histogram index
phot_g_y = []  # Bin index
phot_g_z = []  # Bin values (heights)
phot_bp_x = []
phot_bp_y = []
phot_bp_z = []
phot_rp_x = []
phot_rp_y = []
phot_rp_z = []

for hist_index, hist in enumerate(histograms_g):
    for bin_index, count in enumerate(hist):
        phot_g_x.append(hist_index)  # Histogram number
        phot_g_y.append(bin_index)   # Bin index in histogram
        phot_g_z.append(count)       # Bin value (height)
for hist_index, hist in enumerate(histograms_bp):
    for bin_index, count in enumerate(hist):
        phot_bp_x.append(hist_index)
        phot_bp_y.append(bin_index)
        phot_bp_z.append(count)
for hist_index, hist in enumerate(histograms_rp):
    for bin_index, count in enumerate(hist):
        phot_rp_x.append(hist_index)
        phot_rp_y.append(bin_index)
        phot_rp_z.append(count)

# Create a 3D scatter plot using Plotly
fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=phot_g_x,
    y=phot_g_y,
    z=phot_g_z,
    mode='markers',
    marker=dict(
        size=1,
        color='green',
        opacity=0.6
    ),
))
fig.add_trace(go.Scatter3d(
    x=phot_bp_x,
    y=phot_bp_y,
    z=phot_bp_z,
    mode='markers',
    marker=dict(
        size=1,
        color='lightblue',
        opacity=0.6
    ),
))
fig.add_trace(go.Scatter3d(
    x=phot_rp_x,
    y=phot_rp_y,
    z=phot_rp_z,
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

histogram_num = math.floor(len(histograms_g) * interesting_part)
histogram1_g = histograms_g[histogram_num]
histogram2_g = histograms_g[len(histograms_g)-1 - histogram_num]
bins = np.arange(len(histogram1_g))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_g,
    mode='lines+markers',
    marker=dict(size=4, color='red'),
    line=dict(color='green', width=2),
    name="green left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_g,
    mode='lines+markers',
    marker=dict(size=4, color='lightblue'),
    line=dict(color='green', width=2),
    name="green right"
))

histogram_num = math.floor(len(histograms_bp) * interesting_part)
histogram1_bp = histograms_bp[histogram_num]
histogram2_bp = histograms_bp[len(histograms_bp)-1 - histogram_num]
bins = np.arange(len(histogram1_bp))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_bp,
    mode='lines+markers',
    marker=dict(size=4, color='red'),
    line=dict(color='lightblue', width=2),
    name="blue left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_bp,
    mode='lines+markers',
    marker=dict(size=4, color='green'),
    line=dict(color='lightblue', width=2),
    name="blue right"
))

histogram_num = math.floor(len(histograms_rp) * interesting_part)
histogram1_rp = histograms_rp[histogram_num]
histogram2_rp = histograms_rp[len(histograms_rp)-1 - histogram_num]
bins = np.arange(len(histogram1_rp))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram1_rp,
    mode='lines+markers',
    marker=dict(size=4, color='green'),
    line=dict(color='red', width=2),
    name="red left"
))

fig.add_trace(go.Scattergl(
    x=bins,
    y=histogram2_rp,
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



# on 3rd display it looks like there's no difference, but there is veryyy small difference
# i think i need to calc difference per star, and then display that.
# cuz displaying sum is too tiny to see, the shift is too gentle.
# so i need to use bp_rp, bp_g, g_rp
# starting over...
