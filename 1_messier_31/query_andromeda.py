import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

Gaia.ROW_LIMIT = -1

coord = SkyCoord(ra=10.684, dec=41.269, unit=(u.degree, u.degree), frame='icrs')

# job = Gaia.cone_search_async(coord, radius=u.Quantity(1.5, u.deg), dump_to_file=True, output_format='fits', output_file='andromeda_cone_data.fits')
# results = job.get_results()

results = Gaia.query_object_async(coordinate=coord, width=u.Quantity(3, u.deg), height=u.Quantity(3, u.deg))
results.write('andromeda_square_data.csv', format='csv', overwrite=True)