# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:04:14 2020

@author: Jorge Antonio Matías López
"""

import pandas as pd
import datetime
from pathlib import Path

def lista(ruta = Path.cwd()):
    '''Función que sirve para obtener una lista con los nombres,
       de los archivos contenidos en una carpeta 
    '''
    return [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]

def get_time():
    '''
    Devuelve la fecha y hora actual
    '''
    return datetime.datetime.now().strftime("%A %d %B %Y %I:%M")

def coord_params(X=0, Y=0, x_0=0, y_0=0, dx=1, dy=1):
    '''
    Parametros
    ----------
    X : Tipo: float  
        Coordenada x a transformar. The default is 0.
    Y : Tipo: float 
        Coordenada y a transformar. The default is 0.
    x_0 : Tipo: float   
        Coordenada x de rederencia. The default is 0.
    y_0 : Tipo: float 
        Coordenada y de rederencia. The default is 0.
    dx : Tipo: int 
        Tamaño de la celda en x. The default is 1.
    dy : Tipo: int
        Tamaño de la celda en x. The default is 1.

    Returns
    -------
    row : TYPE int
        Fila correspondiente a la coordenada transformada
    column : TYPE int
        Columna correspondiente a la coordenada transformada
    r_off : TYPE float
    c_off : TYPEfloat
        Offset desde el centro de la celda a la posición real
    '''
    xmod = X - x_0
    ymod = y_0 - Y 
    row = ymod//dy + 1
    column = xmod//dx + 1
    r_off = (ymod - row*dy)/dy
    c_off = (xmod - column*dx)/dx
    return row, column, r_off, c_off

    
def no3_xls2tob(MaxConcObs=200, MaxFluxObs=1, MaxFluxCells=0,#dataset1
                outnam='prueba', inConcObs=21, inFluxObs=22, inSaveObs=23,#dataset2
                Cscale=1.0, iOutCobs=1, iConcLOG=0, iConcINTP=1,#dataset3
                COBSNAM='NO3', layer=4, weight=1, iComp=1,#dataset4
                x0=455204.440, y0=2110063.17+64000, dx=2000, dy=2000, year_ini = 1934,#transformacion coordenadas
                nFluxGroup=1, FScale=1, iOutFlux=1,#dataset6
                nFluxTimeObs=1, ncells=0, iSSType=2,#dataset7
                FluxTimeObs=0, weight_fobs=1, FluxObs=0):#dataset8
    
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
        
        for nombre in nombres:
            df = pd.read_excel(Path.cwd()/ 'input'/'nitratos'/nombre,na_values=0)
            df = df[['ALCALDIA','POZO','X','Y','NO3']]
            cad = ''
            for n in [int(s) for s in nombre if s.isdigit()]:
                cad = cad + str(n)
            df['YR'] = float(cad)
            conjunto_db.append(df)
        dataframe = pd.concat(conjunto_db, ignore_index=True)
        
        nConcObs = dataframe['NO3'].count()
        outnam = outnam + '.tob'
        escritura = open(Path.cwd()/ 'output'/ outnam, 'w') 
        time = get_time()
        header = f'# TOB: Transport Observation package Created on {time} \n'
        escritura.write(header)
        
        data_set1 = f'{MaxConcObs}\t{MaxFluxObs}\t{MaxFluxCells} # Data Set 1: MaxConcObs, MaxFluxObs, MaxFluxCells\n'
        escritura.write(data_set1)
        
        data_set2 = f'{outnam}\t{inConcObs}\t{inFluxObs}\t{inSaveObs} # Data Set 2: OUTNAM, inConcObs, inFluxObs, inSaveObs\n'
        escritura.write(data_set2)
        
        data_set3= f'{nConcObs}\t{Cscale:E}\t{iOutCobs}\t{iConcLOG}\t{iConcINTP}# Data Set 3: nConcObs, CScale, iOutCobs, iConcLOG, iConcINTP\n'
        escritura.write(data_set3)
        check = dataframe['NO3'].notnull()
        
        escritura.write('# Data Set 4: COBSNAM, Layer, Row, Column, iComp, TimeObs, Roff, Coff, weight, COBS\n')
        for fila in range(dataframe.shape[0]):
            if check[fila]:
                COBS = dataframe['NO3'][fila]
                row, column, r_off, c_off = coord_params(dataframe['X'][fila], dataframe['Y'][fila],
                                                         x0, y0, dx, dy)

                #z, rho = geom(row,column) 
                yr = dataframe['YR'][fila]
                timeObs = (yr - year_ini) * 3600 * 24 * 365
                COBS = dataframe['NO3'][fila]# dx*dy*z * rho * (1g/1000mg)
                data_set4 = f'{COBSNAM}\t{layer}\t{int(row)}\t{int(column)}\t{iComp}  {timeObs:E}   {r_off:E}   {c_off:E}   {weight:E}   {COBS:E}\n'#datos del pozo
                escritura.write(data_set4)
        data_set6 = f'\t{nFluxGroup}  {FScale:E}   {iOutFlux} # Data Set 6: nFluxGroup, FScale, iOutFlux\n'
        escritura.write(data_set6)
        data_set7 = f'\t{nFluxTimeObs}    {ncells}    {iSSType} # Data Set 7: nFluxTimeObs, ncells, iSSType\n'
        escritura.write(data_set7)
        FOBSNAM=COBSNAM
        data_set8 = f'{FOBSNAM}   {iComp}   {FluxTimeObs:E}   {weight_fobs:E}   {FluxObs:E} # Data set 8: FOBSNAM, iComp, FluxTimeObs, weight_fobs, FluxObs\n'
        escritura.write(data_set8)
        escritura.close()
                
        return dataframe
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
        encabezado = f'''#Archivo ob_hob\t Created on {get_time()}  
#DATOS DEL ENCABEZADO
#numero_hobs-Número de observaciones de carga 
#ml_obs-Observaciones multicapa  
#max_m-Número máximo de capas para observaciones multicapa
#iu_hobsv-Número de unidad para guardar el archivo de datos de observación
#hob_dry-Valor de el equivalente simulado que es escrito en el archivo de salida
# de observaciones cuando la información es omitida al considerar la celda seca
#tm_of_mult_hbs-Multiplicador de desplazamiento de tiempo para las observaciones de carga(o T/T)\n
'''
        
        formato = f'{numero_hobs} {ml_obs} {max_m} {iu_hobsv} {hob_dry:E}\n{tm_of_mult_hbs:.6f}\n'
        escritura = open(Path.cwd()/ 'output'/'datos.ob_hob', 'w')
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
            row, column, r_off, c_off = coord_params(datos['X'][fila], datos['Y'][fila],
                                                         x0, y0, dx, dy)
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
###########################
datos = no3_xls2tob()
piez = xls2modflow('pozos_excel.xls',455204.44, 2110063.17+64000)


