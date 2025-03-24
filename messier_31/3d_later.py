ra = np.array(results['ra'])
dec = np.array(results['dec'])
distance = 2.537 * 10**6 # TODO for each one

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