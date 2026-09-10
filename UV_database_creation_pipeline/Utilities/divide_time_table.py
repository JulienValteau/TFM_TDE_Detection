import sys
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

def divide_time_table(argv):
    """
    Save the catalogue table time in three subparts
    
    Note that the three value is hard-coded and corresponds to the 
    minimum subparts that has been necessary to have enough RAM  
    to run the match algorithm.

    :arg[1]: Path to the input data to divide 
    :arg[2]: Path to the divided table
    """

    # Open the time table and convert it into pandas
    with fits.open(argv[1]) as hdul:
        df = Table(hdul[-1].data).to_pandas()

    # Split the tables in three parts
    df_chunks = np.array_split(df,3)

    # Save each pandas chunk into single FITS table
    for j, df_chunk in enumerate(df_chunks):
        Table.from_pandas(df_chunk).write(argv[2]+'/'+ 'temp_part_'+str(j+1)+'.fits', overwrite=True)

if __name__ == "__main__":
    divide_time_table(sys.argv)
