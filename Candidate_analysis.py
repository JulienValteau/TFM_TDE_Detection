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
                ind = ind_select2 & (df['CAT'] == cat)
                ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                            yerr=df.loc[ind,flux_err_key], \
                            fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                            ecolor=color_dict[flux])
            ax.errorbar(df.loc[ind,'MJD'],df.loc[ind,flux_key],
                        yerr=df.loc[ind,flux_err_key], \
                        fmt=fmt_dict[df.loc[df.loc[ind,'CAT'].index[0],'CAT']], color=color_dict[flux], \
                        ecolor=color_dict[flux],label=flux)
            
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
        fig.savefig(roots_tfm + TDE_path +'/Test/SCRNUM_CUV_'+str(srcnum)+'.png',dpi=1200)
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
       result1,result2,result3,result4 = run_galex_client_source_tilt(df_results.loc[index,'RA'],
                                                                      df_results.loc[index,'DEC'],
                                                                      60)

       if len(result1)>2 and len(result3)>2 and len(result2)>2 and len(result4)>2:
           if (result1[3].split(' ')[0] == result2[3].split(' ')[0]) and \
           (result2[3].split(' ')[0]==result3[3].split(' ')[0]) and \
           (result3[3].split(' ')[0]==result4[3].split(' ')[0]) and 'No data' not in result1[3]:
               result = result1[3].split(' ')[0]
               print(df_results.loc[index,'SRCNUM_CUV'], ' ',result,' ',Time(result, format="isot", scale="utc").mjd)
               df_cor.loc[index,'MJD'] = Time(result, format="isot", scale="utc").mjd 
           else:
               print('Trying points closer to the source')
               result1,result2,result3,result4 = run_galex_client_source_tilt(df_results.loc[index,'RA'],
                                                                              df_results.loc[index,'DEC'],
                                                                              20)
               if len(result1)>2 and len(result3)>2 and len(result2)>2 and len(result4)>2:
                   if (result1[3].split(' ')[0] == result2[3].split(' ')[0]) and \
                   (result2[3].split(' ')[0]==result3[3].split(' ')[0]) and \
                   (result3[3].split(' ')[0]==result4[3].split(' ')[0] and 'No data' not in result1[3]):
                      
                       result = result1[3].split(' ')[0]
                       print(df_results.loc[index,'SRCNUM_CUV'], ' ',result,' ',Time(result, format="isot", scale="utc").mjd)
                       df_cor.loc[index,'MJD'] = Time(result, format="isot", scale="utc").mjd 
                   else : 
                       if 'No data' in result1[3]:
                           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result1.stdout.split('\n')[3],', ',result2.stdout.split('\n')[3], \
                            ', ',result3.stdout.split('\n')[3],', ' , result4.stdout.split('\n')[3],', ')
                       else:
                           print(df_results.loc[index,'SRCNUM_CUV'],  'NaN : ',result1.stdout.split('\n')[3].split(' ')[0],', ',result2.stdout.split('\n')[3].split(' ')[0], \
                            ', ',result3.stdout.split('\n')[3].split(' ')[0],', ' , result4.stdout.split('\n')[3].split(' ')[0],', ') 
               else:
                   print(result1,result2,result3,result4)
       else:
           print(result1,result2,result3,result4)
    
    # Create the detection flag
    df_cor['DETECTION_FLAG'] = True
     
    # For sources that have no galex point, find upper limit
    ind_galex_source = df_cor['SRCNUM_CUV'].isin(df_test_galex['SRCNUM_CUV'])
    df_no_galex_detection = df_cor.loc[~ind_galex_source,:].groupby('SRCNUM_CUV').agg({'SRCNUM_CUV':'first',\
                                                                                       'RA':'first', \
                                                                                       'DEC' : 'first'})
    for index in df_no_galex_detection.index:
        result= subprocess.run(["python","client_galex.py", 
                            str(round(df_no_galex_detection.loc[index,'RA'],5)), 
                            str(round(df_no_galex_detection.loc[index,'DEC'],5))], 
                           capture_output=True,text=True).stdout.split('\n')
        print(result)
        #df_cor.iloc[-1,['DETECTION_FLAG','CAT','UVW2_FLUX','SRCNUM_CUV','MJD']] = [False, 'GALEX_NUV',results,results]
    
    
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
    result1= subprocess.run(["python","client_galex.py", 
                            str(round(ra+tilt/3600,5)), 
                            str(round(dec,5))], 
                           capture_output=True,text=True).stdout.split('\n')
    result2= subprocess.run(["python","client_galex.py", 
                            str(round(ra-tilt/3600,5)), 
                            str(round(dec,5))], 
                           capture_output=True,text=True).stdout.split('\n')
    result3= subprocess.run(["python","client_galex.py", 
                            str(round(ra,5)), 
                            str(round(dec+tilt/3600,5))], 
                           capture_output=True,text=True).stdout.split('\n')
    result4= subprocess.run(["python","client_galex.py", 
                            str(round(ra,5)), 
                            str(round(dec-tilt/3600,5))], 
                           capture_output=True,text=True).stdout.split('\n')
    
    return result1,result2,result3,result4

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
    TDE_path = '/TDE/Candidates/Ratio higher than 3 and rest lower than 2 (UVW2 and UVM2), UVW1 lower than 2/Best/All bands/Best/Not retained'
    
    ## Print all bands
    list_best_CUV = get_list_best_CUV(TDE_path)
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
