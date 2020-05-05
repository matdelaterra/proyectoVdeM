# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:04:14 2020

@author: Jorge Antonio Matías López
"""
import numpy as np
import pandas as pd
from pathlib import Path

def lista(ruta = Path.cwd()):
    '''Función que sirve para obtener una lista con los nombres
       de los archivos contenidos en una carpeta 
    '''
    return [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]

def traslacion(df=None, x0=0, y0=0, dx=1, dy=1 ):
    ''' Función que sirve para cambiar las coordenadas de UTM a coordenadas locales'''
    if isinstance(df, pd.DataFrame):
        df.loc[:, 'X'] -= x0
        df.loc[:, 'Y'] = y0 - df.loc[:, 'Y']
    
def importar_nitratos():
    
    #seleccionando archivos de excel
    nombres = []
    for archivo in lista(Path.cwd()/ 'input'/'nitratos'):
        if archivo.count('xls') != 0:
            nombres.append(archivo)
        
   
    if len(nombres) != 0:
        
        conjunto_db = []
        #instruccion para leer la geometría
        
        #DataFrame principal
        print('Importando archivos')
        conjunto_db.append(pd.DataFrame({'NOMBRE':pd.Series([]),
                              'CLAVE':pd.Series([]),
                              'X':pd.Series([]), 
                              'Y':pd.Series([]), 
                              'NO3':pd.Series([])}))
        
        for nombre in nombres:
            excel = pd.read_excel(Path.cwd() / 'input' /'nitratos'/ nombre)
            excel = excel[['NOMBRE','CLAVE','X','Y','NO3']]
            conjunto_db.append(excel)
        
        return conjunto_db#pd.concat(conjunto_db)
    else:
        print('No hay archivos para procesar en el directorio')
        

def xls2modflow(archivo=None, x0=0, y0=0, ml_obs=0, max_m=2, iu_hobsv=42, hob_dry=1.0E+30, 
                tm_of_mult_hbs=1.0, na_val = ['ND'], layer=2, year=1934):
    if  isinstance(archivo,str):
        
        datos = pd.read_excel(Path.cwd()/ 'input'/'piezometria'/archivo, na_values=na_val)
        traslacion(datos, x0, y0)
        obs = datos.iloc[:,3:]
        numero_hobs = np.sum(obs.count())        
        formato = f'{numero_hobs} {ml_obs} {max_m} {iu_hobsv} {hob_dry:E}\n{tm_of_mult_hbs:.6f}\n'
    
        escritura = open('datos.ob_hob', 'w')
        escritura.write(formato)
        obs_bool = obs.isnull()
        
        
        for fila in range(obs.shape[0]):
            obsname = datos['ID'][fila]
            x = datos['X'][fila]
            y = datos['Y'][fila]
            hobs = 0.0
            itt = 1
            
            formato = f'{obsname} {layer} {x} {y} {hobs} {itt} \n'
            escritura.write(formato)
            
            for columna in range(obs.shape[1]):
                if not bool(obs_bool.iloc[fila,columna]):
                    yr = list(obs.columns)[columna]
                    obs_subname =  obsname + '_' + str(yr)[-2:]
                    irefsp = 1
                    toffset = (yr - year) * 3600 * 24 * 365
                    hobs = obs.iloc[fila,columna]
                    formato = f'{obs_subname} {irefsp} {toffset:E} {hobs} {itt} \n'
                    escritura.write(formato)
        escritura.close()
    return datos



datos = importar_nitratos()
piez = xls2modflow('pozos_excel.xls',455204.44, 2110063.17+64000)

salida = Path.cwd() / 'output' # ruta a la carpeta de salida para almacenar las gráficas

#obsformation = pd.read_csv(Path.cwd() / 'input' / 'obs_formation.csv') #lectura del csv obs_formation
    
#df = pd.read_fwf(Path.cwd().joinpath('input').joinpath('.hob_out'))




