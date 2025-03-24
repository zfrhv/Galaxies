import numpy as np
from astropy.io import fits
from astropy.table import Table

# import results
hdul = fits.open('andromeda_data.fits')
data = hdul[1].data
results = Table(data)
hdul.close()

useful_columns = []
useless_columns = []
for column in results.columns:
  if np.issubdtype(results[column].dtype, np.floating):
    nan_count = np.sum(np.isnan(results[column]))
  else:
    nan_count = 0
  if nan_count > len(results) * 0.01:
    useless_columns.append(column)
  else:
    useful_columns.append(column)

print('useful columns:')
for column in useful_columns:
  print(column, ":", results[column][0])
