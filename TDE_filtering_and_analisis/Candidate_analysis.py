# This script allows to get deeper analysis of the retained TDE candidates
# In particular, 
#    1- It tries to get the Galex point data if existing or the associated upper-limit
#    2- It computes the W1-W2 from Wise results obtained from Regalade
#    3- It creates the link to go to ESAsky
# All infos are saved in a file TDE_candidates_infos.csv
#
# It prints also the multiband light curves that are saved in the same repository

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gc
#from astroquery.simbad import Simbad
#from astropy import coordinates
#import astropy.units as u
import os
import subprocess
#from matplotlib.widgets import Button
#import webbrowser
#import mpld3
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Angle


def print_light_curve(df,TDE_path):
    """
    This functions allows to print the light curves of various sources

    : df: dataframe of sources observation
    : TDE_path: The path within roots_tfm to save the figures
    """

    fmt_dict = {'GALEX_FUV': 'd', 'GALEX_NUV' : 'o', 'UVOT': '+', 'OM': 'x'}
    color_dict = {'UVW2': 'indigo','UVW1': 'fuchsia','UVM2':'orange','U':'crimson','V':'limegreen','B':'royalblue'}

    # Loop over the unique sources in the provided dataframe
    for srcnum in df['SRCNUM_CUV'].unique():
                       
        fig, ax =plt.subplots()
        ax.set_yscale('log')

        # Store the index of the row from the current source
        ind_select = (df['SRCNUM_CUV']==srcnum)

        # Loop over all bands
        for flux in ['UVW2','UVM2','UVW1','U','V','B']:

            flux_key = flux + "_FLUX"
            flux_err_key = flux_key + "_ERR"
            flux_quality = flux + "_QUALITY_FLAG"

            # From the index keep only the one with the best quality
            ind_select2 = ind_select & (df[flux_quality]==0)

            # Get the list of catalogues
            list_cat = df.loc[ind_select2,'CAT'].unique()           

            # Loop over the catalogues
            for cat in list_cat:

                # If the entries from the catalogue is a detection
                ind = ind_select2 & (df['CAT'] == cat) & (df['DETECTION_FLAG'] == True)
                if len(np.where(ind)[0]) > 0:

                    # Plot it as an errorbar
                    ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                                yerr=df.loc[ind,flux_err_key], \
                                fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                                ecolor=color_dict[flux])
                    
                # If the entries from the catalogue is an upper-limit
                ind2 = ind_select2 & (df['CAT'] == cat) & (df['DETECTION_FLAG'] == False)
                if len(np.where(ind2)[0]) > 0:
                    
                    # Plot it as an errorbar using an upper-limit arrow
                    # TODO: Correct the size of the error so it is fixed and not depends on ymin and ymax
                    ymin, ymax = ax.get_ylim()
                    ax.errorbar(df.loc[ind2,'MJD'],df.loc[ind2,flux_key],
                                yerr=0.05 * (ymax - ymin), \
                                fmt=fmt_dict[df.loc[df.loc[ind2,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                                ecolor=color_dict[flux],uplims=True,capsize=0)

            # Replot the last point to add the label only once
            # TODO: Change the approach to avoid drawing various time the same point potentially with arrows from different size
            if len(np.where(ind)[0]) > 0:
                ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                            yerr=df.loc[ind,flux_err_key], \
                            fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                            ecolor=color_dict[flux],label=flux)
            else:
                ymin, ymax = ax.get_ylim()
                ax.errorbar(df.loc[ind2,'MJD'],df.loc[ind2,flux_key],
                            yerr=0.05 * (ymax - ymin), \
                            fmt=fmt_dict[df.loc[df.loc[ind2,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                            ecolor=color_dict[flux],uplims=True,capsize=0,label=flux)
                
            
        # Add the label, the legend, the text, get the RA, DEC to add to the title and save the figure
        ax.set_xlabel('Time [Modified Julian Date]')
        ax.set_ylabel("Flux [" + "$erg.s^{-1}.cm^{2}.A^{-1}]$")
        ax.legend(loc='best',prop={'size': 8},title= "Flux band")
        ra=df.loc[ind_select.index[np.where(ind_select)][0],'RA']
        dec=df.loc[ind_select.index[np.where(ind_select)][0],'DEC']
        #ax_bouton = plt.axes([0.4, 0.95, 0.2, 0.075])  # position [x, y, largeur, hauteur]
        #bouton = Button(ax_bouton, 'Check on ESASky')
        #bouton.on_clicked(lambda event: open_site(event,ra,dec))
        ax.text(0.6, 0.95, 'GALEX NUV : ' + fmt_dict['GALEX_NUV'] + '   GALEX FUV : ' + fmt_dict['GALEX_FUV'] \
                + '\n     UVOT : ' + fmt_dict['UVOT'] + '         OM : ' + fmt_dict['OM'],\
                transform=ax.transAxes, fontsize=8, verticalalignment='top')
        plt.title("Source (RA DEC) : " + "{0:.5f}".format(ra) + " " + \
                  "{0:.5f}".format(dec))
        fig.subplots_adjust(left=0.25)
        fig.savefig(f"{roots_tfm}{TDE_path}/Test/SCRNUM_CUV_{int(srcnum)}.png",dpi=1200,bbox_inches=None)
        fig.clf()
        plt.close()
        gc.collect()
            
def print_source_infos(df_results,TDE_path):
    """
    This functions allows to print the additional info of the source

    : df_results: dataframe of sources observation
    : TDE_path: The path within roots_tfm to save the csv

    
    : return: df_cor dataframe with additional points (upper-limit), 
    corrected Galex date and a column indicating which point corresponds 
    to a detection.
    """

    # Create the corrected dataframe
    df_cor = df_results.copy()

    # Extracted the source info only from the input dataframe (not the flux)
    df_results_sources = df_results[['SRCNUM_CUV','RA','DEC','logM', 'W1mag','W2mag','D']].drop_duplicates('SRCNUM_CUV').sort_values(['SRCNUM_CUV']).reset_index(drop=True)
    
    ##
    # Addition of galex observation date
    ## 
    
    # Selection of the Galex points from the catalogue
    df_test_galex = df_results[(df_results['CAT']=='GALEX_NUV')]

    # Loop over all the points 
    for index in df_test_galex.index:

        # Get the dates of the observations associated to the 4 points 
        # around the TDE position at 1 arcmin
       result_list = run_galex_client_source_tilt(df_results.loc[index,'RA'],\
                                                  df_results.loc[index,'DEC'],\
                                                  60)
        
       # If at list three points provides info (not NaN)
       if len(np.where([len(result)>2 for result in result_list])[0]) > 2:
           
           # Get the date list, the possible date from the max occurence and the frequence associated
           date_list = [result[3].split(' ')[0] for result in result_list if len(result) > 2]
           possible_date = max(date_list,key=date_list.count)
           frequence = date_list.count(possible_date)

           # If at list three points provide the same date, print it and correct it in the corrected dataframe
           if (frequence>2) and ('No' not in possible_date) and (possible_date):
               print(df_results.loc[index,'SRCNUM_CUV'], ' ',possible_date,' ', 
                     Time(possible_date, format="isot", scale="utc").mjd,' Frequence :',frequence)
               df_cor.loc[index,'MJD'] = Time(possible_date, format="isot", scale="utc").mjd 

           # Otherwise try closer
           else:
               print('Trying points closer to the source')
               
               # Get the dates of the observations associated to the 4 points 
               # around the TDE position at 20 arcsec
               result_list = run_galex_client_source_tilt(df_results.loc[index,'RA'],
                                                          df_results.loc[index,'DEC'],
                                                          20)     

               # If at list three points provides info (not NaN)
               if len(np.where([len(result)>2 for result in result_list])[0]) > 2: 

                   # Get the date list, the possible date from the max occurence and the frequence associated
                   date_list = [result[3].split(' ')[0] for result in result_list if len(result) > 2]
                   possible_date = max(date_list,key=date_list.count)
                   frequence = date_list.count(possible_date)

                   # If at list three points provide the same date, print it and correct it in the corrected dataframe
                   if (frequence>2) and ('No' not in possible_date) and (possible_date): 
                       print(df_results.loc[index,'SRCNUM_CUV'], ' ',possible_date,' ',
                             Time(possible_date, format="isot", scale="utc").mjd,' Frequence :',frequence)
                       df_cor.loc[index,'MJD'] = Time(possible_date, format="isot", scale="utc").mjd 

                   # Otherwise print all info possible 
                   # from the results provided by the client
                   else : 
                        if 'No' in possible_date:
                            print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list[0][3],', ',result_list[1][3], \
                            ', ',result_list[2][3],', ' , result_list[3][3])
                        else:
                           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',date_list[0],', ', 
                                 date_list[1],', ', date_list[2],', ' , 
                                 date_list[3],', ') 
                
               # Otherwise print all info possible 
               # from the results provided by the client
               else:
                   print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list)
           
       # Otherwise print all info possible 
       # from the results provided by the client
       else:
           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list)

    ##
    # Search for upper limit
    ##
    
    # Create the detection flag
    df_cor['DETECTION_FLAG'] = True
     
    # Select source that have no Galex source
    ind_galex_source = df_cor['SRCNUM_CUV'].isin(df_test_galex['SRCNUM_CUV']0)
    df_no_galex_detection = df_cor.loc[~ind_galex_source,:].groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first',\
                                                                                       'RA':'first', \
                                                                                       'DEC' : 'first'})

    # Define the calibration coefficient obtained for Galex against OM bands
    beta_NUV_UVW2 = [-16.51546352,0.34005235,0.10287565]
    beta_NUV_UVM2 =[-16.71595217,0.37372027,0.125861]
    beta_NUV_UVW1 = [-16.95676674,0.36327972,0.2130301]

    # For all index in the dataframe of source with no galex detection
    for index in df_no_galex_detection.index:

        # Run the client rounding the RA,DEC to 6 decimal
        result= subprocess.run(["python","client_galex.py", 
                            str(round(df_no_galex_detection.loc[index,'RA'],6)), 
                            str(round(df_no_galex_detection.loc[index,'DEC'],6))], 
                           capture_output=True,text=True).stdout.split('\n')
        
        # Detect the results containing upper limits 
        if '<' in result[3]:

            # Save the NUV flux in Jansky and print it
            nuv_flux_jy = float(result[3].split(' < ')[1])
            print(df_no_galex_detection.loc[index,'SRCNUM_CUV'],  'Upper limit : ',Time(result[3].split(' < ')[0], format="isot", scale="utc").mjd,' ', 
                  ' NUV:',str(result[3].split(' < ')[1]),
                  ' UVW2:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVW2)),
                  ' UVM2:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVM2)), 
                  ' UVW1:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVW1)))

            # Add the new point to the corrected dataframe, using the calibration to convert the Jansky into erg.s-1.cm2.A-1
            df_cor.loc[len(df_cor),:] = np.nan
            df_cor.loc[len(df_cor)-1,['DETECTION_FLAG','CAT',\
                            'UVW2_FLUX','UVM2_FLUX','UVW1_FLUX',\
                            'UVW2_QUALITY_FLAG','UVM2_QUALITY_FLAG','UVW1_QUALITY_FLAG',\
                            'SRCNUM_CUV','MJD']] = \
                       [False, 'GALEX_NUV',\
                        10**f(np.log10(nuv_flux_jy),beta_NUV_UVW2),\
                        10**f(np.log10(nuv_flux_jy),beta_NUV_UVM2),\
                        10**f(np.log10(nuv_flux_jy),beta_NUV_UVW1),\
                        0, 0, 0, \
                        df_no_galex_detection.loc[index,'SRCNUM_CUV'],\
                        Time(result[3].split(' < ')[0], format="isot", scale="utc").mjd]

        # Otherwise if no data has been detected, print the info
        elif 'No data' in result[3]:
            print(df_no_galex_detection.loc[index,'SRCNUM_CUV'],  'No data : RA ',ra, ' DEC ', dec)
        
        # Otherwise if a point associated to the source have been detected, 
        # compute and print the separation
        else:
            ra = float(result[1].split(':')[1].split(' ')[2])
            dec = float(result[1].split(':')[1].split(' ')[3])
            c1=SkyCoord(ra=df_no_galex_detection.loc[index,'RA']*u.deg, dec=df_no_galex_detection.loc[index,'DEC']*u.deg, frame='icrs')
            c2=SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            print(df_no_galex_detection.loc[index,'SRCNUM_CUV'],  'Detection : RA ',ra, ' DEC ', dec,' Sep ', Angle(str(c1.separation(c2))).arcsec)
    
    
    # Addition of Mass and wise bands
    df_results_sources['diffW1_W2'] = abs(df_results_sources['W1mag']-df_results_sources['W2mag'])
      
    # Addition of ESASky_link column
    for ind in df_results_sources.index:
        df_results_sources.loc[ind,'ESASky_link'] = 'https://sky.esa.int/esasky/?target='+str(df_results_sources.loc[ind,'RA'])+' '+str(df_results_sources.loc[ind,'DEC'])+\
                '&hips=GALEX+GR6%2F7+AIS+color&fov=0.08125700464842057&projection=SIN&cooframe=J2000&sci=true&lang=en'
    
    # Rounding RA and DEC
    df_results_sources['RA'] = round(df_results_sources['RA'],5)
    df_results_sources['DEC'] = round(df_results_sources['DEC'],5)
    
    # Write resulting csv
    df_results_sources.to_csv(roots_tfm + TDE_path + '/TDE_candidates_infos.csv')

    # Return corrected dataframe
    return df_cor
    
        
def get_list_best_CUV(TDE_path):
    """
    This functions allows to read the sources number 
    from the light curves picture saved in a repository

    : TDE_path: The path within roots_tfm where to search for the picture files
    
    : return: list of source number.
    """
                
    list_CUV = []
    local_path = roots_tfm + TDE_path
    contents = os.listdir(local_path)

    # For each element in the repository
    for item in contents:
        item_path = os.path.join(local_path, item)

        # If the element is a file and contains SCRNUM in is name
        if (os.path.isfile(item_path)) and ('SCRNUM' in item) and (len(item.split('_'))>2):

            # Extract the number and add it to the list
            list_CUV.append(np.int32(item.split('_')[2].split('.')[0]))  

    # Return the list of 
    # Combined UV catalogue source number
    return list_CUV

    
def run_galex_client_source_tilt(ra,dec,tilt):
    """
    This functions allows to call the galex client at four points
    around a RA, DEC position

    : ra: Right Ascension in degrees of the consultation
    : dec: Declination in degrees of the consultation
    : tilt: Distance in arcsec with respect to the RA, DEC position to consult
    
    : return: list of the resulting consultations.
    """
    result =[]
    
    # Consultation at RA+tilt/3600
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra+tilt/3600,5)), 
                                  str(round(dec,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))

    # Consultation at RA-tilt/3600
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra-tilt/3600,5)), 
                                  str(round(dec,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))

    # Consultation at DEC+tilt/3600
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra,5)), 
                                  str(round(dec+tilt/3600,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))

    # Consultation at DEC-tilt/3600
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra,5)), 
                                  str(round(dec-tilt/3600,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))

    # Return the list of results of the consultation
    return result


def f(x: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """
    This function defines a second degrees polynom 
    used for the Galex flux calibration

    : x: Flux in input 
    : beta: Calibration coefficient
    
    : return: Calibrated flux
    """
    b0, b1, b2 = beta
    return b0+ b1*x+ b2*x**2

#def open_site(event,ra,dec):
#    webbrowser.open("https://sky.esa.int/esasky/?target="+str(ra)+" "+str(dec)+"&hips=GALEX+GR6%2F7+AIS+color&fov=0.08125700464842057&projection=SIN&cooframe=J2000&sci=true&lang=en")


if __name__ == "__main__":

    """
    This is the main function that perform all the analysis
    """
    
    roots_tfm = '/home/julien/Documents/Etudes/Astrofisica/Master/TFM/Data'
    df_entries = pd.read_csv(roots_tfm + '/Entries_galaxies', low_memory=False)

    # Selected TDE Candidates repository
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2),  UVW1 lower than 2/Best/All bands/Not in other selection'
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2),  UVW1 lower than 2/Best/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 5 and rest lower than 2/Potential TDE/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2), color inversion/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 3 and rest lower than 2 (UVW2 and UVM2), UVW1 lower than 2/Best/All bands/Best/Not retained'
    #TDE_path = '/TDE/Candidates/Ratio higher than 3 and rest lower than 2 (UVW2 and UVM2), UVW1 lower than 2 - rerun/Best/'
    TDE_path = '/TDE/Candidates/Post-filtering-best-candidates'
    #TDE_path = '/TDE/Candidates/Random cases'
    

    # Define the list of source to get the info and print the light curves
    # From the monoband light curves saved in a repository
    #list_best_CUV = get_list_best_CUV(TDE_path)
    # From a manual list
    list_best_CUV = [44098, 45292, 124333, 128873, 190956, 208380, 214117, 229525,
                     230085, 259263]    
    # Random to see how looks random sources in the catalogue
    #list_best_CUV = np.random.randint(1, 500000, 50)

    # Extract the entries associated to the selected sources from the complete database
    df_results = df_entries[df_entries['SRCNUM_CUV'].isin(list_best_CUV)]

    # Removed the Galex FUV points
    df_results = df_results[df_results['CAT'] != 'GALEX_FUV']

    # Get the selected source info, correct the galex date and get the upper limit when available
    df_results_cor_galex = print_source_infos(df_results,TDE_path)

    # Print the all bands light curves from the corrected dataframe
    print_light_curve(df_results_cor_galex,TDE_path)

    # Test to call simbad server
    # df_results_group = df_results.groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first','RA':'first', \
    #                                                          'DEC' : 'first', 'POSERR' : 'first'})
    # simbad = Simbad()
    # simbad.add_votable_fields("otype")
    # df_type = pd.DataFrame(columns=['TYPE','SRCNUM_CUV'])
    # for srcnum in df_results_group['SRCNUM_CUV']:
    #     df_type_src = pd.DataFrame()
    #     ind_select = (df_results_group['SRCNUM_CUV']==srcnum)
    #     c = coordinates.SkyCoord(ra=df_results_group.loc[ind_select,'RA'].values*u.deg, \
    #                              dec=df_results_group.loc[ind_select,'DEC'].values*u.deg, \
    #                                  frame='icrs')
    #     r = 5*u.arcsec
    #     results = simbad.query_region(c, radius=r).to_pandas()["otype"]
    #     if len(results)>0:
    #         df_type_src['TYPE'] = results
    #     else:
    #         df_type_src['TYPE'] = np.nan
    #     df_type_src['SRCNUM_CUV'] = srcnum
    # df_type = pd.concat([df_type, df_type_src], ignore_index=True)
    # df_type.to_csv(roots_tfm + TDE_path +'Candidates_types.csv')
