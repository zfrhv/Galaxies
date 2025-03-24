import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

Gaia.ROW_LIMIT = -1

coord = SkyCoord(ra=10.684, dec=41.269, unit=(u.degree, u.degree), frame='icrs')
job = Gaia.cone_search_async(coord, radius=u.Quantity(1.5, u.deg), dump_to_file=True, output_format='fits', output_file='andromeda_data.fits')
results = job.get_results()