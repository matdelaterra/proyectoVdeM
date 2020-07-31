# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:04:14 2020

@author: Jorge Antonio Matías López
"""
import pandas as pd
import numpy as np
from pathlib import Path
import math
import datetime


def lista(ruta=Path.cwd(), ext='xls' ):
    '''
    Función que sirve para obtener una lista con los nombres 
        de los archivos de determinada extensión contenidos en una carpeta 

    Parámetros
    ----------
    ruta : Path, opcional
        Ruta de la carpeta. Default: Path.cwd().
    ext : str, opcional
        Extensión del archivo. Default: 'xls'.

    Returns
    -------
    nombres : list
        Lista con str de los archivos en el path

    '''
    nombres = []
    for archivo in [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]:
        if archivo.count(ext) != 0:
            nombres.append(archivo)    
    return nombres

def get_time():
    '''
    Devuelve la fecha y hora actual
    '''
    return datetime.datetime.now().strftime("%A %d %B %Y %I:%M")

def coord_params(X=0, Y=0, x_0=0, y_0=0, dx=1, dy=1):
    '''
    Función para transformar las coordenadas utm a coordenadas locales del modelo
    y utilizarlas para obtener la fila, columna y offsets de la observación
    Parametros
    ----------
    X : float  
        Coordenada x a transformar. Default: 0.
    Y : float 
        Coordenada y a transformar. Default: 0.
    x_0 : float   
        Coordenada x de referencia (Esquina superior izquierda). Default: 0.
    y_0 : float 
        Coordenada y de referencia (Esquina superior izquierda). Default: 0.
    dx : int 
        Tamaño de la celda en x. The default is 1.
    dy : int
        Tamaño de la celda en x. The default is 1.

    Returns
    -------
    row : int
        Fila correspondiente a la coordenada transformada
    column : int
        Columna correspondiente a la coordenada transformada
    r_off : float
    c_off : float
        Offset desde el centro de la celda a la posición real
    '''
    xmod = X - x_0
    ymod = y_0 - Y 
    row = int(ymod//dy + 1)
    column = int(xmod//dx + 1)
    r_off = (ymod - row*dy)/dy
    c_off = (xmod - column*dx)/dx
    
    return row, column, r_off, c_off

    
def nitratos2tob(MaxConcObs=1000, MaxFluxObs=1, MaxFluxCells=0,#dataset1
                outnam='prueba', inConcObs=21, inFluxObs=22, inSaveObs=23,#dataset2
                Cscale=1.0, iOutCobs=1, iConcLOG=0, iConcINTP=1,#dataset3
                COBSNAM='NO3', layer=4, weight=1, iComp=1,#dataset4
                x0=455204.440, y0=2110063.17+64000, dx=2000, dy=2000, year_ini = 1934,#transformacion coordenadas
                nFluxGroup=1, FScale=1, iOutFlux=1,#dataset6
                nFluxTimeObs=1, ncells=0, iSSType=2,#dataset7
                FluxTimeObs=0, weight_fobs=1, FluxObs=0,#dataset8
                rho=0.1, na_val = ['ND', 'nd']):
    '''
    Función que convierte una base de datos con información de nitratos 
    de pozos (NO3 en mg/L) a formato .tob (transport observation package)
    de mt3d
    Los parámetros están descritos en la documentación de mt3d
    ---------------------
    Parámetros relevantes
    outnam : str
        Nombre del archivo de salida. Default: 'prueba'
    rho : float
        Porosidad efectiva de la capa seleccionada. Default: 0.1
    x0 : float
        Coordenada x correspondiente al 0 en coordenadas locales. Default: 455204.44 
    y0 : float 
        Coordenada y correspondiente al 0 en coordenadas locales. Default: 2174063.17
    dx : float
        Tamaño de la celda en x. Default: 2000
    dy : float 
        Tamaño de la celda en y. Default: 2000
    year_ini : int
        Año inicial de las observaciones: Default: 1934
    
    '''
    
    #seleccionando archivos de excel
    path_datos = Path.cwd()/ 'tob'/'base de datos'
    path_geom = Path.cwd()/ 'tob'/'geom'
    lista_nombres = lista(path_datos)
    lista_geom = lista(path_geom)

    if len(lista_nombres) != 0 and  len(lista_geom) == 1:
        
        conjunto_db = []
        print('Importando archivos...')
        print('Archivos de base de datos:')

        #DataFrame principal
        for nombre in lista_nombres:
            df = pd.read_excel(path_datos/nombre,na_values=na_val)
            df = df[['ALCALDIA','POZO','X','Y','NO3']]
            fecha = ''
            for n in [int(s) for s in nombre if s.isdigit()]:
                fecha = fecha + str(n)
            df['YR'] = float(fecha)
            conjunto_db.append(df)
            print(nombre)
        dataframe = pd.concat(conjunto_db, ignore_index=True)
        
        #instruccion para leer la geometría 
        geometria = pd.read_excel(path_geom/lista_geom[0])
        print('Archivo de geometría: ', lista_geom[0])
        
        nConcObs = dataframe['NO3'].count()
        nombre_salida = outnam + '.tob'
        escritura = open(Path.cwd()/'tob'/'salida'/nombre_salida, 'w') 
        dtime = get_time()
        header = '# TOB: Transport Observation package\t Created on ' + dtime + ' by Dpto Rec. Nat, IGF UNAM'+ '\n'
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
                yr = dataframe['YR'][fila]
                timeObs = (yr - year_ini) * 3600 * 24 * 365
                
                datos_geom = geometria[(geometria['COLUMN']==column) & (geometria['ROW']==row)]
                dz = datos_geom.iloc[0,4]-datos_geom.iloc[0,5]
                
                COBS = dataframe['NO3'][fila] * dx * dy * dz * rho  # mg/l *(1g/1000mg)(1000l/1m3) *(m3)* (adim) = g
                
                data_set4 = f'{COBSNAM}\t{layer}\t{row}\t{column}\t{iComp}  {timeObs:E}   {r_off:E}   {c_off:E}   {weight:E}   {COBS:E}\n'#datos del pozo
                
                escritura.write(data_set4)               
                
        data_set6 = f'\t{nFluxGroup}  {FScale:E}   {iOutFlux} # Data Set 6: nFluxGroup, FScale, iOutFlux\n'
        escritura.write(data_set6)
        data_set7 = f'\t{nFluxTimeObs}    {ncells}    {iSSType} # Data Set 7: nFluxTimeObs, ncells, iSSType\n'
        escritura.write(data_set7)
        FOBSNAM=COBSNAM
        data_set8 = f'{FOBSNAM}   {iComp}   {FluxTimeObs:E}   {weight_fobs:E}   {FluxObs:E} # Data set 8: FOBSNAM, iComp, FluxTimeObs, weight_fobs, FluxObs\n'
        escritura.write(data_set8)
        escritura.close()
        
        print(f"Finalizado: archivo exportado en {Path.cwd()/'tob'/'salida'/nombre_salida}")
                

    else:
        print('Error en la carga de base de datos y/o geometría')
     



def piezometria2ob_hob(outnam='output', x0=455204.440, y0=2110063.17+64000, dx=2000, dy=2000, ml_obs=0, max_m=2, iu_hobsv=42, hob_dry=1.0E+30, 
                tm_of_mult_hbs=1.0, na_val = ['ND','nd'], layer=2, year=1934):
    ''' 
    Función que convierte un archivo de datos de piezometría al formato ob_hob de modflow
    '''
    
    path_piezometria = Path.cwd()/'ob_hob'/'piezometria'
    lista_piezometria = lista(path_piezometria)
    
    if len(lista_piezometria) == 1:
        #Lectura del excel
        print('Importando archivos...')
        print(f'Cargando datos desde {lista_piezometria[0]}')
        datos = pd.read_excel(path_piezometria/lista_piezometria[0], na_values=na_val)
        obs = datos.iloc[:,3:]
        numero_hobs = obs.count().sum()
        obs_bool = obs.isnull()
        encabezado = f'''#Archivo ob_hob\t Created on {get_time()} by Dpto Rec. Nat, IGF UNAM 
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
        archivo = outnam + '.ob_hob'
        escritura = open(Path.cwd()/ 'ob_hob'/'salida'/archivo, 'w')
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
            obsname = str(datos['ID'][fila])
            row, column, r_off, c_off = coord_params(datos['X'][fila], datos['Y'][fila], x0, y0, dx, dy)
            irefsp = -obs.loc[fila].count()
            toffset = 0.0
            hobs = 0.0
            itt = 1 #1 para usar observaciones de carga, 2 para usar cambios de carga como observaciones
            
            formato = f'{obsname} {layer} {row} {column} {irefsp} {toffset:E} {r_off:E} {c_off:E} {hobs:E}\n{itt}\n'#datos del pozo
            escritura.write(formato)
            
            for columna in range(obs.shape[1]):
                if not obs_bool.iloc[fila,columna]:
                    yr = list(obs.columns)[columna]
                    obs_subname =  obsname + '_' + str(yr)[-2:]
                    irefsp = 1 #periodo de stress para el cual el tiempo se observacion es referenciado
                    toffset = (int(yr) - year) * 3600 * 24 * 365
                    hobs = obs.iloc[fila,columna]
                    formato = f'{obs_subname} {irefsp} {toffset:E} {hobs:E}\n'
                    escritura.write(formato)
        escritura.close()
        print(f"Finalizado: archivo exportado en {Path.cwd()/ 'ob_hob'/'salida'/archivo}")
    else:
        print('Error al cargar el archivo de datos, verifique que se encuentre el archivo correctamente')
        print(f'En el directorio de entrada se encuentran {len(lista_piezometria)} archivos, podría estar abrierto el archivo')



def mk_pfile(outname='pesos', datos=None, sum_sd=None, sum_sd_dr=None):
    '''
    Función de creación de archivos de pesos con formato txt

    Parámetros
    ----------
    outname : str, opcional
        Nombre asignado al archivo de salida. Default: 'pesos'.
    datos : pandas.Dataframe, opcional
        Datos de piezometría. Default: None.
    sum_sd : pandas.Dataframe, opcional
        Datos desviaciones de cargas. Default: None.
    sum_sd_dr : oandas.Dataframe, opcional
        Datos desviaciones de abatimientos. The default is None.

    Returns
    -------
    None.

    '''
    
    if isinstance(datos, pd.DataFrame) and isinstance(sum_sd, pd.DataFrame) and isinstance(sum_sd_dr, pd.DataFrame):
        nombre = outname + '.txt'
        print('Creando archivo: '+nombre)
        archivo = open(Path.cwd()/ 'pesos'/'salida'/nombre, 'w')
        for fila in range(datos.shape[0]):
            pozo = datos['ID'][fila]
            indice = 1
            for columna in range(3,datos.shape[1]):
                if datos.notnull().iloc[fila, columna]:
                    obs_pozo = pozo + '_' + str(indice)
                    if indice == 1:
                        heads = 'Heads'
                        peso = sum_sd.iloc[fila, columna]
                        referencia = datos.iloc[fila, columna]
                        obs = referencia
                    else:
                        heads = 'Head_Changes'
                        peso = sum_sd_dr.iloc[fila, columna]
                        obs = referencia - datos.iloc[fila, columna]

                    archivo.write(f'{obs_pozo} {heads} {obs:.5f} {peso:.5f}\n')
                    indice += 1
        archivo.close()
        print(f"Finalizado: archivo exportado en {Path.cwd()/ 'pesos'/'salida'/nombre}")
        
        
def pesos(nombre='pesos', calcular=False, na_val = ['ND','nd']):
    '''
    Función para crear un fichero con la informacion de los pesos de 
    las observaciones de carga y abatimiento para calibrar en UCODE.
    Permite cargar o calcular los datos a partir de un archivo de excel 
    utilizando la función mk_pfile

    Parámetros
    ----------
    nombre : str, opcional
        Nombre del archivo de salida. Default: 'pesos'.
    calcular : bool, opcional
        Define si se cargarán los datos calculados o se realizarán los cálculos. Default: False.

    Returns
    -------
    None.

    '''
    if calcular:
        path_entrada = Path.cwd()/'pesos'/'datos'/'computar'
        lista_entrada = lista(path_entrada)
        
        if len(lista_entrada) == 1:
            print('Leyendo datos...')
            datos_entrada = pd.ExcelFile(path_entrada/lista_entrada[0])
            hojas = [datos_entrada.parse(hoja, na_values=na_val) for hoja in datos_entrada.sheet_names]
            
            #SD1
            print('Calculando SD1...')
            hojas[1]['sd1'] = np.nan
            for i in range(hojas[1].shape[0]):
                hojas[1].iloc[i, 2] = math.sqrt(((hojas[1].iloc[i, 1]/100*5)/4)**2 + (1/4)**2 + (90*hojas[1].iloc[i, 1]/100/4)**2)
            sd1 = hojas[1].copy()
            sd1.drop('ELEVATION', axis=1, inplace = True)
            sd1.set_index('ID', inplace = True)
            
            #SD2
            print('Calculando SD2...')
            datos = hojas[0].copy()
            grad1 = hojas[2].copy()
            grad2 = hojas[3].copy()

            grad1.set_index('ID', inplace=True)
            grad2.set_index('ID', inplace=True)

            datos['m'] = np.nan
            datos['b'] = np.nan
            for i in range(datos.shape[0]):
                id_pozo =datos.loc[i,'ID']
                datos.loc[i,'m'] = (grad2['ELEVATION'][id_pozo] - grad1['ELEVATION'][id_pozo])/(2005-1984)
                datos.loc[i,'b'] = grad1['ELEVATION'][datos.loc[i,'ID']]-datos.loc[i,'m']*1984
                for j in range(3,datos.shape[1]-2):
                    datos.iloc[i,j] = (datos.loc[i,'m']*datos.columns[j]+datos.loc[i,'b'])/100*90/4
            sd2 = datos.iloc[:,3:datos.shape[1]-2]
            
            
            #sd3
            print('Calculando SD3...')
            df_datos = hojas[0].copy()
            sd3 = df_datos.iloc[:,0:1].copy()
            sd3['sd3'] = np.nan
            niveles = df_datos.iloc[:,3:]
            for fila in range(niveles.shape[0]):
                d = []
                for columna in range(niveles.shape[1]-1):
                   
                    if niveles.notnull().iloc[fila,columna] and niveles.notnull().iloc[fila,columna+1]:
                        d.append(niveles.iloc[fila, columna]-niveles.iloc[fila, columna+1])
            
                if len(d) > 1:
                    d_mean = np.array(d).mean()
                    sd_calc = 0
                    for diff in d:
                        sd_calc += (diff-d_mean)**2
                    sd_calc = math.sqrt(sd_calc/(len(d)-1))
                    sd3.loc[fila,'sd3'] = sd_calc
            sd3.fillna(value=sd3['sd3'].mean())
            sd3.set_index('ID', inplace = True)
            

            #sd4
            print('Calculando SD4...')
            datos1 = hojas[0].copy()
            elevacion = hojas[4].copy()
            elevacion.set_index('ID', inplace=True)

            for i in range(datos.shape[0]):
                datos1.iloc[i,3:] = elevacion['elev pozo'][datos1.loc[i,'ID']] - datos1.iloc[i,3:]

            sd4 = datos1.iloc[:,3:]
            sd4 = (sd4 * 0.001) / 2
            sd4.fillna(value=0, inplace=True)
            sd4
            
            #sd5 
            print('Calculando SD5...')
            datos.iloc[:,:-2] = hojas[0]
            for i in range(datos.shape[0]):
                for j in range(3,datos.shape[1]-2):
                    datos.iloc[i,j] = (datos.loc[i,'m']*datos.columns[j]+datos.loc[i,'b'])/100*2000/4
            sd5 = datos.iloc[:,3:datos.shape[1]-2]
            
            #SdCarga
            print('Calculando suma sd carga...')
            sum_sd_carga = hojas[0].copy()
            for i in range(sum_sd_carga.shape[0]):
                for j in range(3,sum_sd_carga.shape[1]):
                    pozo = sum_sd_carga.loc[i,'ID']
                    sum_sd_carga.iloc[i,j] = math.sqrt(sd1['sd1'][pozo]**2 + sd2.iloc[i, j-3]**2 + sd3['sd3'][pozo]**2
                                                       + sd4.iloc[i, j-3]**2 + sd5.iloc[i, j-3]**2)
                    
            #Suma Abatimiento
            print('Calculando suma sd abatimiento...')
            sum_sd_abat = hojas[0].copy()
            for i in range(sum_sd_abat.shape[0]):
                for j in range(3,sum_sd_abat.shape[1]):
                    pozo = sum_sd_abat.loc[i,'ID']
                    sum_sd_abat.iloc[i,j] = math.sqrt(sd3['sd3'][pozo]**2 + sd4.iloc[i, j-3]**2)
            
            mk_pfile(nombre,hojas[0], sum_sd_carga,sum_sd_abat)
        else:
            print(f'Error al cargar archivo de entrada, verificar que se encuentre en {path_entrada}')
    else:
        path_entrada = Path.cwd()/'pesos'/'datos'/'completo'
        lista_entrada = lista(path_entrada)
        if len(lista_entrada) == 1:
            print('Importando archivos...')
            print(f'Cargando datos desde {lista_entrada[0]}')
            excel = pd.ExcelFile(path_entrada/lista_entrada[0])
            datos = excel.parse(excel.sheet_names[0],na_values=na_val)
            suma_sd = excel.parse(excel.sheet_names[-2])
            suma_sd_dr = excel.parse(excel.sheet_names[-1])
            mk_pfile(nombre, datos, suma_sd, suma_sd_dr)


    

