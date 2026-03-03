import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
roots_tfm = '/home/julien/Documents/Etudes/Astrofisica/Master/TFM/Data'

def divide_time_table(argv):
    """
    Save the catalogue table time in various subpart

    :arg[0]: Repository in Data/ of the catalogue 
    :arg[1]: Name of the catalogue file
    """
    with fits.open(roots_tfm + '/'+argv[1]+'/'+argv[2]) as hdul:
        df = Table(hdul[2].data).to_pandas()

    df_chunks = np.array_split(df,3)

    for j, df_chunk in enumerate(df_chunks):
        Table.from_pandas(df_chunk).write(roots_tfm + '/'+argv[1]+'/'+ 'temp_part_'+str(j+1)+'.fits', overwrite=True)

if __name__ == "__main__":
    divide_time_table(sys.argv)
