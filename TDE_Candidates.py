import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gc
#from astroquery.simbad import Simbad
#from astropy import coordinates
#import astropy.units as u

def select_candidates_ratio(df_entries,bands,threshold,test_type):
    list_flux_dict = {}
    for flux in bands:
        flux_name = flux + '_FLUX'
        flux_err = flux + '_FLUX_ERR'
        flux_quality = flux + '_QUALITY_FLAG'
        
        df_entries_flux = df_entries[['SRCNUM_CUV','MJD',flux_name,flux_err,flux_quality,'CAT']].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_entries_flux = df_entries_flux[(df_entries_flux[flux_quality]==0) & (df_entries_flux['CAT']!='GALEX_FUV')]
        
        df_entries_flux['COUNT'] = df_entries_flux['SRCNUM_CUV'].copy()
        df_entries_flux_group = df_entries_flux.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})
        
        if test_type == 'equal':
            df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>3]['SRCNUM_CUV'])]
        elif test_type == 'superior':
            df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>1]['SRCNUM_CUV'])]
        else:
            raise ValueError(str(test_type)+'is not a specific test type. Choose between ''equal''or ''superior''')
        
        df_entries_flux[flux_name +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_name].shift(-1)
        df_entries_flux[flux_err +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_err].shift(-1)
        df_entries_flux['MJD_next'] = df_entries_flux.groupby('SRCNUM_CUV')['MJD'].shift(-1)

        df_entries_flux[flux_name +'_var'] = df_entries_flux[flux_name +'_next'] / df_entries_flux[flux_name]

        df_entries_test = df_entries_flux[(df_entries_flux[flux_name +'_var']>threshold)]

        #print(df_entries_flux[[flux_name,'MJD',flux_name +'_var']])  
        df_entries_test['COUNT'] = df_entries_test['SRCNUM_CUV'].copy()
        df_entries_test_groups=df_entries_test.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count',\
                                                                          flux_name +'_var': 'min',flux_name +'_var': 'max'})
        if test_type == 'equal':
            ind_select = (df_entries_test_groups['COUNT']==1)
        elif test_type == 'superior':
            ind_select = (df_entries_test_groups['COUNT']>=1)

        list_flux_dict[flux_name] = pd.DataFrame() 
        list_flux_dict[flux_name]=df_entries_test_groups.loc[ind_select,:] 

        print("Number of candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))

    return list_flux_dict


def select_candidates_diff(df_entries):
    list_flux_dict = {}
    for flux in ['UVW2','UVM2','UVW1']:
        flux_name = flux + '_FLUX'
        flux_err = flux + '_FLUX_ERR'
        flux_quality = flux + '_QUALITY_FLAG'
        
        df_entries_flux = df_entries[['SRCNUM_CUV','MJD',flux_name,flux_err,flux_quality,'CAT']].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_entries_flux = df_entries_flux[(df_entries_flux[flux_quality]==0) & (df_entries_flux['CAT']!='GALEX_FUV')]
        
        df_entries_flux['COUNT'] = df_entries_flux['SRCNUM_CUV'].copy()
        df_entries_flux_group = df_entries_flux.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count'})
        df_entries_flux = df_entries_flux[df_entries_flux['SRCNUM_CUV'].isin(df_entries_flux_group[df_entries_flux_group['COUNT']>2]['SRCNUM_CUV'])]

        
        df_entries_flux[flux_name +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_name].shift(-1)
        df_entries_flux[flux_err +'_next'] = df_entries_flux.groupby('SRCNUM_CUV')[flux_err].shift(-1)
        df_entries_flux['MJD_next'] = df_entries_flux.groupby('SRCNUM_CUV')['MJD'].shift(-1)

        df_entries_flux[flux_name +'_var'] = (df_entries_flux[flux_name] - df_entries_flux[flux_name +'_next']) / \
        np.sqrt( df_entries_flux[flux_err]**2 + df_entries_flux[flux_err +'_next']**2)

        df_entries_test = df_entries_flux[(df_entries_flux[flux_name +'_var']<-30)]

        df_entries_test['COUNT'] = df_entries_test['SRCNUM_CUV'].copy()
        df_entries_test_groups=df_entries_test.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first', 'COUNT' : 'count',\
                                                                      flux_name +'_var': 'min',flux_name +'_var': 'max'})
          
        ind_select = (df_entries_test_groups['COUNT']==1)

        list_flux_dict[flux_name] = pd.DataFrame() 
        list_flux_dict[flux_name]=df_entries_test_groups.loc[ind_select,:] 

        print("Number of candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))

    return list_flux_dict


def print_light_curve_by_flux(df,flux):
    fmt_dict = {'GALEX_FUV': 'd', 'GALEX_NUV' : 'o', 'UVOT': '+', 'OM': 'x'}
    color_dict = {'UVW2': 'indigo','UVW1': 'darkviolet','UVM2':'deeppink','U':'crimson','V':'limegreen','B':'royalblue'}
        
    for srcnum in df['SRCNUM_CUV'].unique():
               
        flux_key = flux + "_FLUX"
        flux_err_key = flux_key + "_ERR"
        flux_quality = flux + "_QUALITY_FLAG"
        
        ind_select = (df['SRCNUM_CUV']==srcnum) & (df[flux_quality]==0)
        list_cat = df.loc[ind_select,'CAT'].unique()

        if len(list_cat)>0:
            
            fig, ax =plt.subplots()
            ax.set_yscale('log')
            
            for cat in list_cat:
                ind = ind_select & (df['CAT'] == cat)
                ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                            yerr=df.loc[ind,flux_err_key], \
                            fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                            ecolor=color_dict[flux],label=cat)

            #print(ind_select)
            for i in ind_select.index[np.where(ind_select)]:
                #print(i)
                ax.annotate(f'{df.loc[i,'OBSID']}',xy=(df.loc[i,'MJD'],df.loc[i,flux_key]),xytext=(0,4),textcoords='offset points',
                            ha='center', va='bottom', size=6)
            
            #ax.axvline(x=df_tde.loc[tde,'MJD'],color='orange',label='TDE detection date')
            ax.set_xlabel('Time [Modified Julian Date]')
            ax.set_ylabel("Flux [" + "$erg.s^{-1}.cm^{2}.A^{-1}]$")
            ax.legend(loc='lower right',prop={'size': 8},title= "Flux band " +flux)
            plt.title("Source (RA DEC) : " + "{0:.5f}".format(df.loc[ind_select.index[np.where(ind_select)][0],'RA']) + " " + \
                      "{0:.5f}".format(df.loc[ind_select.index[np.where(ind_select)][0],'DEC']))
            fig.savefig(roots_tfm + '/TDE/Candidates/SCRNUM_CUV_'+str(srcnum)+'_'+flux+'_flux.png',dpi=1200)
            fig.clf()
            plt.close()
            gc.collect()

if __name__ == "__main__":

    roots_tfm = '/home/julien/Documents/Etudes/Astrofisica/Master/TFM/Data'
    df_entries = pd.read_csv(roots_tfm + '/Entries_galaxies_more_than_2.csv')
    
    #df_entries = df_entries[df_entries['SRCNUM_CUV']==1544]
    
    #list_flux_dict = select_candidates_diff(df_entries)
    
    # Search for candidates where one point has a variance ratio higher than 4 and all 
    # the other values have variance ratio lower than 1.5 in band UVW2 and UVM2
    bands = ['UVW2', 'UVM2']
    list_flux_dict_up=select_candidates_ratio(df_entries,bands,4,'equal')    
    list_flux_dict_down=select_candidates_ratio(df_entries,bands,2,'equal')
    
    # Search for NO candidates where at least one point has a variance ratio higher than 2 in 
    # UVW1 band
    bands = ['UVW1']
    list_flux_dict_no_candidates=select_candidates_ratio(df_entries,bands,2,'superior')
    
    list_flux_dict={}
    for flux in ['UVW2','UVM2']:
        flux_name = flux+'_FLUX'
        list_flux_dict[flux_name] = pd.DataFrame() 
        list_flux_dict[flux_name]['SRCNUM_CUV'] = list_flux_dict_up[flux_name].loc[list_flux_dict_up[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_down[flux_name]['SRCNUM_CUV']),'SRCNUM_CUV']
        list_flux_dict[flux_name]['SRCNUM_CUV'] = list_flux_dict[flux_name].loc[~list_flux_dict[flux_name]['SRCNUM_CUV'].isin(list_flux_dict_no_candidates['UVW1_FLUX']['SRCNUM_CUV']),'SRCNUM_CUV']
        print("Number of final candidates in band " + flux_name) 
        print(str(list_flux_dict[flux_name].shape[0]))
  
    for flux in ['UVW2','UVM2']:
        flux_name = flux+'_FLUX'
        flux_quality = flux + '_QUALITY_FLAG'
        df_last_entries = pd.read_csv(roots_tfm + '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2),  UVW1 lower than 2/Candidates_'+flux+'.csv')
        df_entries2= df_entries[~df_entries['SRCNUM_CUV'].isin(df_last_entries['SRCNUM_CUV'])]
        df_results = df_entries2[df_entries2['SRCNUM_CUV'].isin(list_flux_dict[flux_name]['SRCNUM_CUV'])].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        #df_results = df_entries[df_entries['SRCNUM_CUV'].isin(list_flux_dict[flux_name]['SRCNUM_CUV'])].dropna(subset=[flux_name]).sort_values(['SRCNUM_CUV','MJD']).reset_index(drop=True)
        df_results = df_results[(df_results[flux_quality]==0) & (df_results['CAT']!='GALEX_FUV')]
        df_results.to_csv(roots_tfm + '/TDE/Candidates/Candidates_'+flux+'.csv')
        print_light_curve_by_flux(df_results,flux)
        
        #df_results_source = df_results[(df_results['SRCNUM_CUV']==1544) | (df_results['SRCNUM_CUV']==202300)]
        #print_light_curve_by_flux(df_results_source,flux)
    
    # df_results_group = df_results.groupby('SRCNUM_CUV').agg({'RA':'first', 'DEC' : 'first', 'POSERR' : 'first'})
    # simbad = Simbad()
    # simbad.add_votable_fields("otype")
    # c = coordinates.SkyCoord(ra=list(df_results_group['RA'])*u.deg, dec=list(df_results_group['DEC'])*u.deg, frame='icrs')
    # r = list(df_results_group['POSERR']*10)* u.arcsec
    # results = simbad.query_region(c, radius=r).to_pandas()
