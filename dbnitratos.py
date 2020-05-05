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
        

def xls2modflow(archivo=None, x0=0, y0=0, dx=2000, dy=2000, ml_obs=0, max_m=2, iu_hobsv=42, hob_dry=1.0E+30, 
                tm_of_mult_hbs=1.0, na_val = ['ND'], layer=2, year=1934):
    
    
    if  isinstance(archivo,str):
        #Lectura del excel
        datos = pd.read_excel(Path.cwd()/ 'input'/'piezometria'/archivo, na_values=na_val)
        obs = datos.iloc[:,3:]
        numero_hobs = obs.count().sum()
        obs_bool = obs.isnull()
        encabezado = '''#DATOS DEL ENCABEZADO
#numero_hobs-Número de observaciones de carga 
#ml_obs-Observaciones multicapa  
#max_m-Número máximo de capas para observaciones multicapa
#iu_hobsv-Número de unidad para guardar el archivo de datos de observación
#hob_dry-Valor de el equivalente simulado que es escrito en el archivo de salida
# de observaciones cuando la información es omitida al considerar la celda seca
#tm_of_mult_hbs-Multiplicador de desplazamiento de tiempo para las observaciones de carga(o T/T)\n
'''
        
        formato = f'{numero_hobs} {ml_obs} {max_m} {iu_hobsv} {hob_dry:E}\n{tm_of_mult_hbs:.6f}\n'
        escritura = open('datos.ob_hob', 'w')
        escritura.write(encabezado)
        encabezado = '''#DATOS DE OBSERVACION
#ENCABEZADO
#obsname-Nombre del pozo
#layer-Capa
#row-Fila
#column-Columna
#irefsp-periodos de stress cuyo tiempo de observacion es referenciado
#toffset-time_offset
#r_offset
#c_offset
#hobs-observacion(0.0)
#itt-1 para usar observaciones de carga, 2 para usar cambios de carga como observaciones
#OBSERVACION
#obs_subname-Nombre de observación
#irefsp-Periodo de stress
#toffset-time offset
#hobs-Observación\n
'''
        escritura.write(encabezado)   
        escritura.write(formato)
        
        for fila in range(obs.shape[0]):
            obsname = datos['ID'][fila]
            xmod = datos['X'][fila] - x0
            ymod = y0 - datos['Y'][fila]
            row = ymod//dy + 1
            column = xmod//dx + 1
            r_off = (ymod - row*dy)/dy
            c_off = (xmod - column*dx )/dx
            
            irefsp = -obs.loc[fila].count()
            toffset = 0.0
            hobs = 0.0
            itt = 1 #1 para usar observaciones de carga, 2 para usar cambios de carga como observaciones
            
            formato = f'{obsname} {layer} {row} {column} {irefsp} {toffset:E} {r_off:E} {c_off:E} {hobs:E}\n{itt}\n'#datos del pozo
            escritura.write(formato)
            
            for columna in range(obs.shape[1]):
                if not bool(obs_bool.iloc[fila,columna]):
                    yr = list(obs.columns)[columna]
                    obs_subname =  obsname + '_' + str(yr)[-2:]
                    irefsp = 1 #periodo de stress para el cual el tiempo se observacion es referenciado
                    toffset = (yr - year) * 3600 * 24 * 365
                    hobs = obs.iloc[fila,columna]
                    formato = f'{obs_subname} {irefsp} {toffset:E} {hobs:E}\n'
                    escritura.write(formato)
        escritura.close()
    return datos



datos = importar_nitratos()
piez = xls2modflow('pozos_excel.xls',455204.44, 2110063.17+64000)




