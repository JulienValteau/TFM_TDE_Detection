import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import pandas

class OneCatalog_Source:
    #Class corresponding to sources from individual catalogs
    def __init__(self, info_dict):
        #Method to build the catalogs sources. When reading the catalog files,
        #the data is stored in a dictionary info_dict, which is read by this
        # function to intialise the object. We can add new info to this dic
        self.src_id = info_dict['src_id']
        self.ra = info_dict['ra']
        self.dec = info_dict['dec']
        self.tab_fluxes = info_dict['tab_fluxes']
        self.tab_dates = info_dict['tab_dates']
        # For upper limits, what to add ?

class MultiCatalog_Source:
    #Class corresponding to sources from merged catalogs - ideally, correspond to one single astrophysical object
    def __init__(self, info_dict):
        self.merged_ra = info_dict['merged_ra']
        self.merged_dec = info_dict['merged_dec']
        self.tab_sources = info_dict['tab_sources']
        self.variability = np.nan

    def plot_lightcurves(self):
        #Complete this function to plot the fluxes of all OneCatalog_Source associated with this MultiCatalog_Source
        #You can loop to access the sources: for instance, "for catalog_src in self.tab_sources:"
        plt.figure()
        for i,catalog_src in enumerate(self.tab_sources):
          plt.scatter(catalog_src.tab_dates,catalog_src.tab_fluxes,color='C'+str(i),label='source #' + str(i))

def read_OM(filename):
    om_catalog = []
    #For OM: the doc is here https://www.cosmos.esa.int/web/xmm-newton/om-catalogue
    # and the data here https://nxsa.esac.esa.int/catalogues/XMM-OM-SUSS6.2.fits.tar.gz
    # Once you downloaded it, you can upload it to the Colab using the folder icon on the left
    with fits.open(filename, memmap=True) as data:
        df = pandas.DataFrame(data[1].data)

        #You will need to
        # 1. Split the dataframe by SRCNUM, so you have all detections of one single source grouped together
        # 2. Modify the class OneCatalog_Source to now load in an object, on top of what we already have:
            # - the source IAUName
            # - the different fluxes in the 6 bands instead of just 1 number
            # - the errors on these fluxes, and on the RA and Dec as well
        # 3. Load each splitted dataframe as a OneCatalog_Source object, and add it to the om_catalog table
    print(len(df))
    return om_catalog
        plt.legend()
