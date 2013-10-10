
import kplr
import pyfits
import numpy as np

client = kplr.API()
dir(client)

# Kepler ID
#  (if of interest) : KOI number
#  (if confirmed)   : "Name" e.g., '62b'

#koi = client.koi(952.01)
planet = client.planet("62b")  # or "Kepler-62 b"

# This object has a lot of attributes (with names given
# http://archive.stsci.edu/search_fields.php?mission=kepler_candidates
# such as a period:

print(planet.koi_period)
# 5.715

# For some reason, the KOI table tends to have more precise measurements so
# we can look at that instead:

koi = planet.koi
print("{0.koi_period} ± {0.koi_period_err1}".format(koi))
# 5.71493 ± 0.00019

# The attributes of the KOI object are given in the `MAST description of the
# kepler/koi table
# <http://archive.stsci.edu/search_fields.php?mission=kepler_koi>`_.
# You can also directly query the KOI table using:

koi = client.koi("256.01")

# To find all light curves
# short_cadence = False (no short cadence data)
# fetch = True (download now, else they are downloaded when first opened)
lightcurves = koi.get_light_curves(short_cadence=False, fetch=True)
for lc in lightcurves:
    print(lc.filename)

# This will download the FITS files containing the light curves to the directory
# given by the ``KPLR_DATA_DIR`` environment variable (or ``~/.kplr/data`` by
# default). To load one of the files, you'll need to make sure that you have
# `pyfits <http://pythonhosted.org/pyfits/>`_ installed and then you can use the
# ``Dataset`` object:

# extract data:
time, flux, ferr, quality = np.array([]), np.array([]), np.array([]), np.array([])
for lc in lightcurves:
    with lc.open(clobber=True) as tmpfile:
        # The lightcurve data are in the first FITS HDU.
        hdu_data = tmpfile[1].data
        time = np.append(time,hdu_data["time"])
        flux = np.append(flux,hdu_data["sap_flux"])
        ferr = np.append(ferr,hdu_data["sap_flux_err"])
        quality = np.append(quality,hdu_data["sap_quality"])
        print "Finished extracting data for file: " + str(lc.filename)


plot(flux)


