#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 22:33:31 2026

@author: julien
"""

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

    fmt_dict = {'GALEX_FUV': 'd', 'GALEX_NUV' : 'o', 'UVOT': '+', 'OM': 'x'}
    color_dict = {'UVW2': 'indigo','UVW1': 'fuchsia','UVM2':'orange','U':'crimson','V':'limegreen','B':'royalblue'}
    #color_dict = {'UVW2': 'C0','UVW1': 'C1','UVM2':'C2','U':'C3','V':'limegreen','B':'royalblue'}
        
    for srcnum in df['SRCNUM_CUV'].unique():
                       
        fig, ax =plt.subplots()
        ax.set_yscale('log')
        ind_select = (df['SRCNUM_CUV']==srcnum)

        for flux in ['UVW2','UVM2','UVW1','U','V','B']:

            flux_key = flux + "_FLUX"
            flux_err_key = flux_key + "_ERR"
            flux_quality = flux + "_QUALITY_FLAG"
            
            ind_select2 = ind_select & (df[flux_quality]==0)
            list_cat = df.loc[ind_select2,'CAT'].unique()           
            
            for cat in list_cat:
                ind = ind_select2 & (df['CAT'] == cat) & (df['DETECTION_FLAG'] == True)
                if len(np.where(ind)[0]) > 0:
                    ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                                yerr=df.loc[ind,flux_err_key], \
                                fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                                ecolor=color_dict[flux])
                ind2 = ind_select2 & (df['CAT'] == cat) & (df['DETECTION_FLAG'] == False)
                if len(np.where(ind2)[0]) > 0:
                    ymin, ymax = ax.get_ylim()
                    ax.errorbar(df.loc[ind2,'MJD'],df.loc[ind2,flux_key],
                                yerr=0.05 * (ymax - ymin), \
                                fmt=fmt_dict[df.loc[df.loc[ind2,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                                ecolor=color_dict[flux],uplims=True,capsize=0)
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
                
            
        #ax.axvline(x=df_tde.loc[tde,'MJD'],color='orange',label='TDE detection date')
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
        plt.tight_layout()
        fig.savefig(f"{roots_tfm}{TDE_path}/Test/SCRNUM_CUV_{int(srcnum)}.png",dpi=1200)
        #mpld3.save_html(fig, roots_tfm + TDE_path +'SCRNUM_CUV_'+str(srcnum)+'.html')
        fig.clf()
        plt.close()
        gc.collect()
            
def print_source_infos(df_results,TDE_path):
    
    df_cor = df_results.copy()
    df_results_sources = df_results[['SRCNUM_CUV','RA','DEC','logM', 'W1mag','W2mag']].drop_duplicates('SRCNUM_CUV').sort_values(['SRCNUM_CUV']).reset_index(drop=True)
    
    # Addition of galex observation date
    df_test_galex = df_results[(df_results['CAT']=='GALEX_NUV')]
    for index in df_test_galex.index:
       result_list = run_galex_client_source_tilt(df_results.loc[index,'RA'],\
                                                  df_results.loc[index,'DEC'],\
                                                  60)
       if len(np.where([len(result)>2 for result in result_list])[0]) > 2:
           date_list = [result[3].split(' ')[0] for result in result_list if len(result) > 2]
           possible_date = max(date_list,key=date_list.count)
           frequence = date_list.count(possible_date)
             
           if frequence>2 and 'No' not in possible_date:
               print(df_results.loc[index,'SRCNUM_CUV'], ' ',possible_date,' ', 
                     Time(possible_date, format="isot", scale="utc").mjd,' Frequence :',frequence)
               df_cor.loc[index,'MJD'] = Time(possible_date, format="isot", scale="utc").mjd 
           else:
               print('Trying points closer to the source')
               result_list = run_galex_client_source_tilt(df_results.loc[index,'RA'],
                                                          df_results.loc[index,'DEC'],
                                                          20)         
               if len(np.where([len(result)>2 for result in result_list])[0]) > 2:                   
                   date_list = [result[3].split(' ')[0] for result in result_list if len(result) > 2]
                   possible_date = max(date_list,key=date_list.count)
                   frequence = date_list.count(possible_date)
                   
                   if frequence>2 and 'No' not in possible_date: 
                       print(df_results.loc[index,'SRCNUM_CUV'], ' ',possible_date,' ',
                             Time(possible_date, format="isot", scale="utc").mjd,' Frequence :',frequence)
                       df_cor.loc[index,'MJD'] = Time(possible_date, format="isot", scale="utc").mjd 
                   else : 
                        if 'No' in possible_date:
                            print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list[0][3],', ',result_list[1][3], \
                            ', ',result_list[2][3],', ' , result_list[3][3])
                        else:
                           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',date_list[0],', ', 
                                 date_list[1],', ', date_list[2],', ' , 
                                 date_list[3],', ') 
               else:
                   print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list)
       else:
           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result_list)
    
    # Create the detection flag
    df_cor['DETECTION_FLAG'] = True
     
    # For sources that have no galex point, find upper limit
    ind_galex_source = df_cor['SRCNUM_CUV'].isin(df_test_galex['SRCNUM_CUV'])
    df_no_galex_detection = df_cor.loc[~ind_galex_source,:].groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first',\
                                                                                       'RA':'first', \
                                                                                      'DEC' : 'first'})
    beta_NUV_UVW2 = [-16.51546352,0.34005235,0.10287565]
    beta_NUV_UVM2 =[-16.71595217,0.37372027,0.125861]
    beta_NUV_UVW1 = [-16.95676674,0.36327972,0.2130301]
    for index in df_no_galex_detection.index:
        result= subprocess.run(["python","client_galex.py", 
                            str(round(df_no_galex_detection.loc[index,'RA'],6)), 
                            str(round(df_no_galex_detection.loc[index,'DEC'],6))], 
                           capture_output=True,text=True).stdout.split('\n')
        #print(result)
        if '<' in result[3]:
            nuv_flux_jy = float(result[3].split(' < ')[1])
            print(df_no_galex_detection.loc[index,'SRCNUM_CUV'],  'Upper limit : ',Time(result[3].split(' < ')[0], format="isot", scale="utc").mjd,' ', 
                  ' NUV:',str(result[3].split(' < ')[1]),
                  ' UVW2:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVW2)),
                  ' UVM2:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVM2)), 
                  ' UVW1:',str(10**f(np.log10(nuv_flux_jy),beta_NUV_UVW1)))
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
        else:
            ra = float(result[1].split(':')[1].split(' ')[2])
            dec = float(result[1].split(':')[1].split(' ')[3])
            c1=SkyCoord(ra=df_no_galex_detection.loc[index,'RA']*u.deg, dec=df_no_galex_detection.loc[index,'DEC']*u.deg, frame='icrs')
            c2=SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            print(df_no_galex_detection.loc[index,'SRCNUM_CUV'],  'Detection : RA ',ra, ' DEC ', dec,' Sep ', Angle(str(c1.separation(c2))).arcsec)
    
    
    # Addition of Mass and wise bands
    df_results_sources['diffW1_W2'] = abs(df_results_sources['W1mag']-df_results_sources['W2mag'])
#     df_results_sources[['LogMass', 'W1','W2']] = np.nan
#     for flux in ['UVW2','UVM2']:
# #        subprocess.run(["bash", "get_regalade_data.sh", flux, str(TDE_path)], 
# #                       capture_output=False,text=True)
#         df_match = pd.read_csv(roots_tfm + TDE_path + '/All bands/Candidates_regalade_info.csv')
#         ind_select = df_results_sources['SRCNUM_CUV'].isin(df_match['SRCNUM_CUV'])
#         df_results_sources.loc[ind_select,'LogMass', 'W1','W2'] = df_match['LogMass', 'W1','W2']
      
    # Addition of ESASky_link
    for ind in df_results_sources.index:
        df_results_sources.loc[ind,'ESASky_link'] = 'https://sky.esa.int/esasky/?target='+str(df_results_sources.loc[ind,'RA'])+' '+str(df_results_sources.loc[ind,'DEC'])+\
                '&hips=GALEX+GR6%2F7+AIS+color&fov=0.08125700464842057&projection=SIN&cooframe=J2000&sci=true&lang=en'
    
    # Rounding RA and DEC
    df_results_sources['RA'] = round(df_results_sources['RA'],5)
    df_results_sources['DEC'] = round(df_results_sources['DEC'],5)
    
    # Write results
    df_results_sources.to_csv(roots_tfm + TDE_path + '/TDE_candidates_infos.csv')
    
    return df_cor
    
        
def get_list_best_CUV(TDE_path):
                
    list_CUV = []
    local_path = roots_tfm + TDE_path
    contents = os.listdir(local_path)
    for item in contents:
        item_path = os.path.join(local_path, item)
        if (os.path.isfile(item_path)) and ('SCRNUM' in item) and (len(item.split('_'))>2):
            list_CUV.append(np.int32(item.split('_')[2].split('.')[0]))  
    
    return list_CUV

def run_galex_client_source_tilt(ra,dec,tilt):
    result =[]
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra+tilt/3600,5)), 
                                  str(round(dec,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra-tilt/3600,5)), 
                                  str(round(dec,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra,5)), 
                                  str(round(dec+tilt/3600,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))
    result.append(subprocess.run(["python","client_galex.py", 
                                  str(round(ra,5)), 
                                  str(round(dec-tilt/3600,5))], 
                                 capture_output=True,text=True).stdout.split('\n'))
    
    return result

def f(x: np.ndarray, beta: np.ndarray) -> np.ndarray:
    b0, b1, b2 = beta
    return b0+ b1*x+ b2*x**2

# Fonction appelée lors du clic
#def open_site(event,ra,dec):
#    webbrowser.open("https://sky.esa.int/esasky/?target="+str(ra)+" "+str(dec)+"&hips=GALEX+GR6%2F7+AIS+color&fov=0.08125700464842057&projection=SIN&cooframe=J2000&sci=true&lang=en")
        
if __name__ == "__main__":

    roots_tfm = '/home/julien/Documents/Etudes/Astrofisica/Master/TFM/Data'
    df_entries = pd.read_csv(roots_tfm + '/Entries_galaxies_more_than_2.csv')
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2),  UVW1 lower than 2/Best/All bands/Not in other selection'
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2),  UVW1 lower than 2/Best/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 5 and rest lower than 2/Potential TDE/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 4 and rest lower than 2 (UVW2 and UVM2), color inversion/All bands/Best'
    #TDE_path = '/TDE/Candidates/Ratio higher than 3 and rest lower than 2 (UVW2 and UVM2), UVW1 lower than 2/Best/All bands/Best/Not retained'
    TDE_path = '/TDE/Candidates/Reanalize_galex'
    
    ## Print all bands
    #list_best_CUV = get_list_best_CUV(TDE_path)
    #list_best_CUV = [8555, 96685, 99220, 154467, 172161, 200050, 227141, 259263,
    #                 276861, 34929, 44098, 124333, 190956, 214117]
    list_best_CUV = [151515]
    df_results = df_entries[df_entries['SRCNUM_CUV'].isin(list_best_CUV)]
    df_results = df_results[df_results['CAT'] != 'GALEX_FUV']
    
    df_results_cor_galex = print_source_infos(df_results,TDE_path)
    
    print_light_curve(df_results_cor_galex,TDE_path)

    
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
