# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:24:13 2019

@author: Jorge Antonio Matías López
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pickle
from pathlib import Path

def graficasN_D(s_name, zone=True, obsinfo_loaded=True, united_graph = True, timestep='s', startdate='1935-01-01', ddn_lim=[-50, 20]):
    '''
    Permite graficar niveles y descensos de un modelo de Modflow a partir de los archivos de salida .hob_out, 
    y un objeto pickle creado a partir del archivo .ob_hob, con la posibilidad de obtener gráficas 
    que muestren los niveles y descensos observados y simulados en un único gráfico,
    o de manera separada obteniendo únicamente gráficos de los datos observados
    Parámetros
    • s_name: Nombre del modelo utilizado(str)
    • zone: Asigna zonas con el archivo obs_formation.csv (opcional, bool, default=True)
    • obsinfo_loaded: Indica si la información hob del modelo ya ha sido obtenida en el objeto pickle(opcional, bool, default=True)
    • united_graph: Define si los resultados serán gráficas unidas o separadas (opcional, bool, default= True)
    • timestep: Unidad de tiempo del modelo (str, default ='s'(segundos))
    • startdate: Fecha inicial del modelo (str, 'aaaa-mm-dd', default= '1935-01-01')
    • ddn_lim: Límites superior e inferior para el gráfico de abatimiento (list,default=[-50, 20])
    '''
    filelocation = Path.cwd() / 'output' # ruta a la carpeta de salida para almacenar las gáficas
    
    '''
    Importa observaciones de carga desde el archivo .hob.out el cual proveé valores simulados y observados de carga
    '''
    if zone:
        geology = ['Lacustrine','Alluvial','Basalt','Volcaniclastic','Andesite']#lista de geología
        obsformation = pd.read_csv(Path.cwd() / 'input' / 'obs_formation.csv') #lectura del csv obs_formation
    
    df = pd.read_fwf(Path.cwd().joinpath('input').joinpath(s_name+'.hob_out'),widths=[22,19,22])#lectura de el archivo hob
    df.columns = ['simulated','observed','obs_name']#cambio de nombre de las columnas
    # Se añaden nuevas columnas con valores NaN
    df['time_series'] = np.nan
    df['obs_id'] = np.nan
    
    if zone:
        df['dem'] = np.nan
        df['geo'] = np.nan
    df['abssimulated'] = np.nan
    df['absobserved'] = np.nan
    
    #se cargan los datos de observación
    if obsinfo_loaded:
        print('Abriendo fichero de información...')
        with open(Path.cwd() / 'input' / 'OBS.pickle', 'rb') as handle:
            obsinfo = pickle.load(handle)
    else:
        print('Procesando archivo de observaviones...')
        mf = flopy.modflow.Modflow.load(Path.cwd().joinpath('input').joinpath(s_name+'.nam'), verbose=True,
                               version='mf2005')
        hob = flopy.modflow.ModflowHob.load(Path.cwd().joinpath('input').joinpath(s_name+'.ob_hob'), mf)
        winfofile = Path.cwd() / 'input' / 'OBS.pickle'
        with open(winfofile, 'wb') as handle:
            pickle.dump(hob.obs_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        obsinfo = hob.obs_data
                
           
    for i in obsinfo:# se recorre cada objeto y se obtiene su información
        
        oname = i.obsname# se guarda el id en oname
        t = np.ones(len(df[df['obs_name'].str.contains(oname)].index))*np.nan 
        t[:i.nobs] = [x[0] for x in i.time_series_data]
        df.loc[df['obs_name'].str.contains(oname),'time_series'] = t  
        df.loc[df['obs_name'].str.contains(oname),'obs_id'] = oname 
        df.loc[df['obs_name'].str.contains(oname),'abssimulated'] = df[df['obs_name'].str.contains(oname)][1:]['simulated'] + df[df['obs_name'].str.contains(oname)]['simulated'].values[0]
        df.loc[df['obs_name']==(oname),'abssimulated'] = df[df['obs_name']==oname]['simulated']##??? 
        df.loc[df['obs_name'].str.contains(oname),'absobserved'] = df[df['obs_name'].str.contains(oname)][1:]['observed'] + df[df['obs_name'].str.contains(oname)]['observed'].values[0]
        df.loc[df['obs_name']==(oname),'absobserved'] = df[df['obs_name']==oname]['observed']#???

    #ciclo para guardar los datos de niveles que son tomados como referencia
    for i in range(df.shape[0]):
        if df['observed'][i] > 1000:
            df.loc[i,'absobserved'] = df['observed'][i]
            df.loc[i, 'abssimulated'] = df['simulated'][i]
    
    if zone:        
        #ciclo para asignar la elevación y el tipo de geología
        for i, r in obsformation.iterrows():
            
            df.loc[df['obs_id']==r['IDPOZO'],'dem'] = r['ELEV']
            df.loc[df['obs_id']==r['IDPOZO'],'geo'] = geology[r['GEOLOGY_ZONE']-1] 

    '''
    Hasta aquí se ha preparado el DataFrame para realizar los recortes para las graficas
    
    '''     
    #Para graficar se utiliza una lista de todos los pozos que se encuentran el en DataFrame
    #En caso de querer graficar algunos pozos en específico solo es necesario modificar esta lista
    l1 = list(df['obs_id'])
    pozos = list(set(list(df['obs_id'])))
    
    #### Registro de los pozos que tienen un dato de niveles y los que tienen descensos también
    un_dato = []
    completo = []
    
    #una opción para probar el funcionamiento del código con unas cuantas gráfica 
    #Para graficar todos comentar esta variable y la condición al final del ciclo
    #e = 0 
    print('Generando gráficas..')
    for i in pozos:
        n = l1.count(i)
        ddn_data = df[df['obs_id']==i].copy()
        ddn_data['time_series'] = pd.to_timedelta(ddn_data['time_series'], timestep) + pd.to_datetime(startdate)        
        ddn_data = ddn_data.set_index('time_series')
        if n == 1:### Únicamente niveles
            
            un_dato.append(i)
    
            if united_graph:
                ### Gráficas unidas
                fig, axes = plt.subplots(figsize=(5,6))
                ddn_data['abssimulated'].plot(ax=axes, color="red", lw=1, ls='-', marker='s',  markersize=4)
                ddn_data['absobserved'].plot(ax=axes, color="blue", lw=1, ls='-', marker='o', markersize=4)
                axes.legend(['Simulado','Observado'])
                axes.set_title('Pozo '+ i)
                axes.xaxis.label.set_visible(False)               
                
                plt.tight_layout()
                filename = filelocation.joinpath(i+'.png')
                #print(filename)
                plt.savefig(str(filename), dpi=300)
                plt.close()
                
                
                    
            else:
                #opcion gráficos separados
                fig1, axes1 = plt.subplots(figsize=(5,6))
                ddn_data['absobserved'].plot(ax=axes1, color="blue", lw=1, ls='-', marker='o', markersize=4)
                axes1.legend(['Observado']) # Leyendas
                axes1.set_title('Niveles Pozo'+ i) # Título
                axes1.xaxis.label.set_visible(False)
                
                plt.tight_layout()
                filename = filelocation.joinpath(i+'.png')
                #print(filename)
                plt.savefig(str(filename), dpi=300)
                plt.close()
                
                     
        else:## Niveles y descensos ######
            
            completo.append(i)
            ddn_data['simulated'].iloc[0] = 0
            ddn_data['observed'].iloc[0] = 0
            
            if united_graph:
                ### Gráficas unidas
                fig, axes = plt.subplots(1, 2, figsize=(12,5))
                ddn_data['abssimulated'].plot(ax=axes[0], color="red", lw=1, ls='-', marker='s',  markersize=4)
                ddn_data['absobserved'].plot(ax=axes[0], color="blue", lw=1, ls='-', marker='o', markersize=4)
                axes[0].legend(['Simulados','Observados'])
                axes[0].set_title('Niveles Pozo '+ i)
                axes[0].xaxis.label.set_visible(False)
                    
                ddn_data['simulated'].plot(ax=axes[1], color="red", lw=1, ls='-', marker='s', markersize=4)
                ddn_data['observed'].plot(ax=axes[1], color="blue", lw=1, ls='-', marker='o', markersize=4)
                axes[1].set_ylim(ddn_lim)
                axes[1].set_title('Descensos Pozo '+ i)
                axes[1].legend(['Simulados','Observados'])
                axes[1].xaxis.label.set_visible(False)
               
                
                plt.tight_layout()
                filename = filelocation.joinpath(i+'.png')
                #print(filename)
                plt.savefig(str(filename), dpi=300)
                plt.close()

                
            else:##############################
                #opcion gráficos separados
                fig1, axes1 = plt.subplots(figsize=(5,6))
                
                ddn_data['absobserved'].plot(ax=axes1, color="blue", lw=1, ls='-', marker='o', markersize=4)
                #ddn_data['abssimulated'].plot(ax=axes1, color="blue", lw=1, ls='--', marker='^')
                axes1.legend(['Observados'])
                axes1.set_title('Niveles Pozo '+ i)## Título
                axes1.xaxis.label.set_visible(False)
                
                fig, axes = plt.subplots(figsize=(5,6))
                ddn_data['observed'].plot(ax=axes,color="blue", lw=1, ls='-', marker='^', markersize=4)
                axes.set_ylim(ddn_lim)
                axes.legend(['Observado'])
                axes.set_title('Descensos Pozo '+ i)#Título
                axes.xaxis.label.set_visible(False)
            
            
                filename_niv = filelocation.joinpath(i+'_niv'+'.png')
                filename_des = filelocation.joinpath(i+'_des'+'.png')
                fig1.savefig(str(filename_niv), dpi=300)
                fig.savefig(str(filename_des), dpi=300)
                plt.close()
                plt.close()
            #e += 1
    
            #if e>4:
             #   break
    print('Finalizado')
    
