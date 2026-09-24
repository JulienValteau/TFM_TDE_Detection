# This script allows to filter the TDE candidates within the UV galaxies
# database loaded as a dataframe (see function main at the end)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gc
#from astroquery.simbad import Simbad
#from astropy import coordinates
#import astropy.units as u

def select_candidates_color(df_entries,bands):
    """
    This functions allows to select candidates according the UV color change

    : df_entries: dataframe over which the filtering is performed
    : bands: list of bands to consider to perform the filtering
    
    """
    
    list_flux_dict = {}
    for flux in bands:
        flux_name = flux + '_FLUX'
        flux_err = flux + '_FLUX_ERR'
        flux_quality = flux + '_QUALITY_FLAG'

        # Create the list of color to consider for the filtering
        # Note that only the difference with respect the less energetic 
        # is considered
        if flux == 'UVW2':
            list_color = [flux+'-UVW1', flux+'-UVM2']
        elif flux =='UVM2':
            list_color = [flux+'-UVW1']
        else:
            raise ValueError(str(flux)+'is not a valid band for this selection. Choose between ''UVW2'' or ''UVM2''')

        list_flux_dict[flux_name] = pd.DataFrame() 

        # For each UV color (UVW2-UVM2 or UVW2-UVW1 or UVM2-UVW1)
        for color in list_color:
            
            df_color = pd.DataFrame()
            
            # Select the entries columns, eliminating the NaN rows and sorting the values according the source number and the observation date
            # Only the best quality entries are kept and the NUV values for Galex.
            df_entries_flux = df_entries[['SRCNUM_CUV','MJD',flux_name,flux_err,flux_quality,color,'CAT']].dropna(subset=[flux_name, color]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
            df_entries_flux = df_entries_flux[(df_entries_flux[flux_quality]==0) & (df_entries_flux['CAT']!='GALEX_FUV')]

            # Create a count column and groupby in order to filters only the sources with more than 2 entries
            df_entries_flux['COUNT'] = df_entries_flux['SRCNUM_CUV'].copy()
            df_entries_flux_group = df_entries_flux.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})
            df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>2]['SRCNUM_CUV'])]

            # Create a column of the next values, shifting the color rows and 
            # compute the variability as the product of the current color and the next one
            df_entries_flux[color +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[color].shift(-1)
            df_entries_flux[color +'_var'] = df_entries_flux[color +'_next'] * df_entries_flux[color]

            # Detect the rows were the sign of the color change
            df_entries_test = df_entries_flux[(df_entries_flux[color +'_var']<0)]

            # Select the sources where the color sign change only once
            df_entries_test['COUNT'] = df_entries_test['SRCNUM_CUV'].copy()
            df_entries_test_groups=df_entries_test.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})            
            ind_select = (df_entries_test_groups['COUNT']==1)
            df_color = df_entries_test_groups.loc[ind_select,:] 

            # Concate to the list already existing for this flux the list of sources selected
            list_flux_dict[flux_name]=pd.concat([list_flux_dict[flux_name],df_color]) 

        print("Number of candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))
    
    return list_flux_dict

def select_candidates_ratio(df_entries,bands,threshold,test_type):
    """
    This functions allows to select candidates according the change 
    of their flux ratio

    : df_entries: dataframe over which the filtering is performed
    : bands: list of bands to consider to perform the filtering
    : threshold: limit for the ratio variation detection
    : test_type: string stating the kind of filtering desired (equal, 'superior or inferior)
    
    """
    list_flux_dict = {}
    for flux in bands:
        flux_name = flux + '_FLUX'
        flux_err = flux + '_FLUX_ERR'
        flux_quality = flux + '_QUALITY_FLAG'

        # Select the entries columns, eliminating the NaN rows and sorting the values according the source number and the observation date
        # Only the best quality entries are kept and the NUV values for Galex.
        df_entries_flux = df_entries[['SRCNUM_CUV','MJD',flux_name,flux_err,flux_quality,'CAT']].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_entries_flux = df_entries_flux[(df_entries_flux[flux_quality]==0) & (df_entries_flux['CAT']!='GALEX_FUV')]

        # Create a count column and groupby in order to filters only the sources with more points than the value fixed
        df_entries_flux['COUNT'] = df_entries_flux['SRCNUM_CUV'].copy()
        df_entries_flux_group = df_entries_flux.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})

        # If we want to select the sources that exceed the fixed variability threshold exactly once, it will need 
        # to have more than three points to detect a peak of variability and be relevant as TDE candidate (4 points => 3 variabilities, a peak
        # is made at least of this combination: low variability point - high variability point - low variability point)
        if test_type == 'equal':
            df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>3]['SRCNUM_CUV'])]
        # If we want to select the sources that comply with the threshold (upper or lower bound), even one variability is enough (at least
        # two points)
        elif (test_type == 'superior') or (test_type == 'inferior') :
            df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>1]['SRCNUM_CUV'])]
        else:
            raise ValueError(str(test_type)+'is not a specific test type. Choose between ''equal'',  ''superior'' or ''inferior''')

        # Create a column of the next values for the flux and compute the ratio dividing
        # the next values by the previous one, so that the peak is positive for TDE
        df_entries_flux[flux_name +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_name].shift(-1)
        df_entries_flux[flux_name +'_var'] = df_entries_flux[flux_name +'_next'] / df_entries_flux[flux_name]

        # Select the sources that exceed the threshold (either upper or lower bound)
        # Note that the searched peaks of variability for TDE are positive, so that 
        # the equal options implies to exceed an upper bound
        if (test_type == 'equal') or (test_type == 'superior'):
            df_entries_test = df_entries_flux[(df_entries_flux[flux_name +'_var']>threshold)]
        elif (test_type == 'inferior'):
            df_entries_test = df_entries_flux[(df_entries_flux[flux_name +'_var']<threshold)]

        # From the sources that exceed the threshold, select the one that exceed it only once for the equal option
        # and at least once for the superior and inferior options
        df_entries_test['COUNT'] = df_entries_test['SRCNUM_CUV'].copy()
        df_entries_test_groups=df_entries_test.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count',\
                                                                          flux_name +'_var': 'min',flux_name +'_var': 'max'})
        if test_type == 'equal':
            ind_select = (df_entries_test_groups['COUNT']==1)
        elif (test_type == 'superior') or (test_type == 'inferior') :
            ind_select = (df_entries_test_groups['COUNT']>=1)

        # Assign to this flux the selected sources dataframe 
        list_flux_dict[flux_name] = pd.DataFrame() 
        list_flux_dict[flux_name]=df_entries_test_groups.loc[ind_select,:] 

        print("Number of candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))

    return list_flux_dict


def select_candidates_diff(df_entries):
    """
    This functions allows to select candidates according the change 
    of their absolute flux

    : df_entries: dataframe over which the filtering is performed    
    """
    list_flux_dict = {}
    for flux in ['UVW2','UVM2','UVW1']:
        flux_name = flux + '_FLUX'
        flux_err = flux + '_FLUX_ERR'
        flux_quality = flux + '_QUALITY_FLAG'

        # Select the entries columns, eliminating the NaN rows and sorting the values according the source number and the observation date
        # Only the best quality entries are kept and the NUV values for Galex.
        df_entries_flux = df_entries[['SRCNUM_CUV','MJD',flux_name,flux_err,flux_quality,'CAT']].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_entries_flux = df_entries_flux[(df_entries_flux[flux_quality]==0) & (df_entries_flux['CAT']!='GALEX_FUV')]

        # Create a count column and groupby in order to filters only for the sources with more than 2 entries
        df_entries_flux['COUNT'] = df_entries_flux['SRCNUM_CUV'].copy()
        df_entries_flux_group = df_entries_flux.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})
        df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>2]['SRCNUM_CUV'])]

        # Create a column of the next values for the flux and the error, shifting the row from one index and  
        # compute the variability as the difference of their flux weighted by the square root of the sum of their quadratic errors
        df_entries_flux[flux_name +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_name].shift(-1)
        df_entries_flux[flux_err +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_err].shift(-1)
        df_entries_flux[flux_name +'_var'] = (df_entries_flux[flux_name] - df_entries_flux[flux_name +'_next']) / \
        np.sqrt( df_entries_flux[flux_err]**2 + df_entries_flux[flux_err +'_next']**2)

        # Select the sources where the variability is lower than -30 (holistic value from known TDE analysis)
        # Note that positive peaks correspond to negative variability)
        df_entries_test = df_entries_flux[(df_entries_flux[flux_name +'_var']<-30)]

        # From the sources where the variability is lower than -30, select the ones that have only one positive peak
        df_entries_test['COUNT'] = df_entries_test['SRCNUM_CUV'].copy()
        df_entries_test_groups=df_entries_test.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count',\
                                                                      flux_name +'_var': 'min',flux_name +'_var': 'max'})      
        ind_select = (df_entries_test_groups['COUNT']==1)

        # Assign to this flux the selected sources dataframe 
        list_flux_dict[flux_name] = pd.DataFrame() 
        list_flux_dict[flux_name]=df_entries_test_groups.loc[ind_select,:] 

        print("Number of candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))

    return list_flux_dict


def print_light_curve_by_flux(df,flux):
    """
    This functions allows to print the light curves of various sources
    for a single flux

    : df: dataframe of sources observation
    : flux: selected flux between UVW2, UVM2, UVW1, U, V, B
    """
    fmt_dict = {'GALEX_FUV': 'd', 'GALEX_NUV' : 'o', 'UVOT': '+', 'OM': 'x'}
    color_dict = {'UVW2': 'indigo','UVW1': 'darkviolet','UVM2':'deeppink','U':'crimson','V':'limegreen','B':'royalblue'}

    # Loop over the unique sources in the provided dataframe
    for srcnum in df['SRCNUM_CUV'].unique():
               
        flux_key = flux + "_FLUX"
        flux_err_key = flux_key + "_ERR"
        flux_quality = flux + "_QUALITY_FLAG"

        # From the current source select only the best quality observations
        ind_select = (df['SRCNUM_CUV']==srcnum) & (df[flux_quality]==0)

        # Get the list of catalogues associated to the current source
        list_cat = df.loc[ind_select,'CAT'].unique()

        # If the list of catalogues is not empty
        # Note: Protection in case no observation from quality 0 is present in the dataframe
        if len(list_cat)>0:
            
            fig, ax =plt.subplots()
            ax.set_yscale('log')

            # Loop over the catalogues
            for cat in list_cat:

                # Draw the point from the selected catalogue
                ind = ind_select & (df['CAT'] == cat)
                ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                            yerr=df.loc[ind,flux_err_key], \
                            fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                            ecolor=color_dict[flux],label=cat)

            # Add the OBSID to each drawn points
            for i in ind_select.index[np.where(ind_select)]:
                ax.annotate(f'{df.loc[i,'OBSID']}',xy=(df.loc[i,'MJD'],df.loc[i,flux_key]),xytext=(0,4),textcoords='offset points',
                            ha='center', va='bottom', size=6)

            # Add the label, the legend, the title, save the figure and clean the memory garbage
            ax.set_xlabel('Time [Modified Julian Date]')
            ax.set_ylabel("Flux [" + "$erg.s^{-1}.cm^{2}.A^{-1}]$")
            ax.legend(loc='lower right',prop={'size': 8},title= "Flux band " +flux)
            plt.title("Source (RA DEC) : " + "{0:.5f}".format(df.loc[ind_select.index[np.where(ind_select)][0],'RA']) + " " + \
                      "{0:.5f}".format(df.loc[ind_select.index[np.where(ind_select)][0],'DEC']))
            #plt.show()
            fig.savefig(roots_tfm + '/TDE/Candidates/SCRNUM_CUV_'+str(srcnum)+'_'+flux+'_flux.png',dpi=1200)
            fig.clf()
            plt.close()
            gc.collect()

if __name__ == "__main__":
    """
    This is the main function that allows to perform the filtering and draw the 
    filtered sources corresponding single flux light curve in order to allow visual inspection
    """

    roots_tfm = '/home/julien/Documents/Etudes/Astrofisica/Master/TFM/Data'
    df_entries = pd.read_csv(roots_tfm + '/Entries_galaxies.csv')

    ## 1- Examples of filtering

    # Difference filtering
    #list_flux_dict = select_candidates_diff(df_entries)

    # Ratio filtering
    # Search for candidates where one point has a variance ratio higher than 3 and all 
    # the other values have variance ratio lower than 2 in band UVW2 and UVM2
    # (Search for candidates with exactly one point with variance than 3 and one
    # point higher than 2 => It will be the same point, all other points have 
    # variance lower than 2 and the point with variance higher, will have a 
    # variance higher than 3)
    bands = ['UVW2', 'UVM2']
    list_flux_dict_up=select_candidates_ratio(df_entries,bands,3,'equal')    
    list_flux_dict_down=select_candidates_ratio(df_entries,bands,2,'equal')

    # Color filtering
    # Search for candidates where the color change of sign
    #df_entries['UVW2-UVW1'] = df_entries['UVW2_FLUX'] - df_entries['UVW1_FLUX']
    #df_entries['UVM2-UVW1'] = df_entries['UVM2_FLUX'] - df_entries['UVW1_FLUX']
    #df_entries['UVW2-UVM2'] = df_entries['UVW2_FLUX'] - df_entries['UVW1_FLUX']
    #list_flux_dict_change = select_candidates_color(df_entries,bands)
    
    # Search for NO candidates where at least one point has a variance ratio higher than 2 in 
    # UVW1 band
    bands = ['UVW1']
    list_flux_dict_no_candidates=select_candidates_ratio(df_entries,bands,2,'superior')
   
    # Note : This is an improvement to eliminate candidate not stable in UVW1 with negative peak
    # Search for NO candidates where at least one point has a variance ratio lower than 0.2 (negative peak) in 
    # all bands
    #bands = ['UVW1']
    #list_flux_dict_no_candidates2=select_candidates_ratio(df_entries,bands,0.2,'inferior')

    # Combine the various list to obtained the desired candidates
    list_flux_dict={}
    
    # Search for candidates selecting according the UVW2 or UVM2 peak
    for flux in ['UVW2','UVM2']:
        flux_name = flux+'_FLUX'
        list_flux_dict[flux_name] = pd.DataFrame() 

        # 1- List of candidates sources number that has at the same time only one point that exceed the first threshold and 
        # one point that exceed a second lower threshold (so that all other point have variability under this second threshold)
        list_flux_dict[flux_name]['SRCNUM_CUV'] = list_flux_dict_up[flux_name].loc[list_flux_dict_up[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_down[flux_name]['SRCNUM_CUV']),'SRCNUM_CUV']

        # 2 - List of candidates sources number that are not in the list of no candidates according the criteria applied to UVW1
        list_flux_dict[flux_name]['SRCNUM_CUV'] = list_flux_dict[flux_name].loc[~list_flux_dict[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_no_candidates['UVW1_FLUX']['SRCNUM_CUV']),'SRCNUM_CUV']

        # (Improvement) 2 - List of candidates sources number that are not in the list of no candidates according the criteria applied to UVW1
        #list_flux_dict[flux_name]['SRCNUM_CUV'] = list_flux_dict[flux_name].loc[(~list_flux_dict[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_no_candidates['UVW1_FLUX']['SRCNUM_CUV'])) & (~list_flux_dict[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_no_candidates2['UVW1_FLUX']['SRCNUM_CUV'])),'SRCNUM_CUV']
      
        print("Number of final candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))

    ## 2- Light curves draw according the selection performed
    for flux in ['UVW2','UVM2']:
        flux_name = flux+'_FLUX'
        flux_quality = flux + '_QUALITY_FLAG'

        # Load the dataframe obtained from the last filtering in order to discard the common candidates
        #df_last_entries = pd.read_csv(roots_tfm + '/TDE/Candidates/Ratio higher than 3 and rest lower than 2 (UVW2 and UVM2), UVW1 lower than 2/Candidates_'+flux+'.csv')
        #df_entries2= df_entries[~df_entries['SRCNUM_CUV'].isin(df_last_entries['SRCNUM_CUV'])]
        #df_results = df_entries2[df_entries2['SRCNUM_CUV'].isin(list_flux_dict[flux_name]['SRCNUM_CUV'])].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)

        # Select the from the entries the one that are in the final selected list of candidates eliminating the NaN rows and 
        # sorting the values according the source number and the observation date
        # Only the best quality entries are kept and the NUV values for Galex.
        df_results = df_entries[df_entries['SRCNUM_CUV'].isin(list_flux_dict[flux_name]['SRCNUM_CUV'])].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_results = df_results[(df_results[flux_quality]==0) & (df_results['CAT']!='GALEX_FUV')]

        # Save the final list in a csv to be able to reload it if necessary
        df_results.to_csv(roots_tfm + '/TDE/Candidates/Candidates_'+flux+'.csv')

        # Print the light curves for the final list
        print_light_curve_by_flux(df_results,flux)
