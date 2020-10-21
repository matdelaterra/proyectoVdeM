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


def corregir_nombre(cadena):
    '''
    Función que recibe una cadena de caracteres, verifica si el primer elemento
    es un carácter alfabético y genera una nueva cadena.
    Los elementos no alfanuméricos son sustituídos por '_'

    Parámetros
    ----------
    cadena : str
        Cadena de caracteres a convertir

    Returns
    -------
    cadena : str
        Cadena original
    n_nombre : TYPE
        Nueva cadena

    '''
    cadena = str(cadena)
    lista = []
    
    if not cadena[0].isalpha():
        lista.append('p')

    for letra in cadena:    
        if letra.isalnum():
            lista.append(letra)

        else:
            lista.append('_')
    n_nombre = ''.join(lista)
    return cadena, n_nombre


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

    
def nitratos2tob(MaxConcObs=1000, MaxFluxObs=0, MaxFluxCells=0,#dataset1
                outnam='prueba', inConcObs=29, inFluxObs=0, inSaveObs=89,#dataset2
                Cscale=1.0, iOutCobs=1, iConcLOG=0, iConcINTP=1,#dataset3
                layer=4, weight=1, iComp=1,#dataset4
                x0=455204.440, y0=2174063.17, dx=2000, dy=2000, year_ini = 1934,#transformacion coordenadas
                nFluxGroup=0, FScale=1, iOutFlux=1,#dataset6
                nFluxTimeObs=0, ncells=0, iSSType=2,#dataset7
                FluxTimeObs=0, weight_fobs=1, FluxObs=0,#dataset8
                rho=0.3, na_val = ['ND', 'nd']):
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
        datos = pd.read_excel(path_datos/lista_nombres[0])## Leer
        datos = datos.iloc[:, 2:]
         
        #instruccion para leer la geometría 
        geometria = pd.read_excel(path_geom/lista_geom[0])
        print('Archivo de geometría: ', lista_geom[0])
        
        nombre_salida = outnam + '.tob'
        escritura = open(Path.cwd()/'tob'/'salida'/nombre_salida, 'w') 
        dtime = get_time()
        header = '# TOB: Transport Observation package\t Created on ' + dtime + ' by Dpto Rec. Nat, IGF UNAM'+ '\n'
        escritura.write(header)
        
        MaxConcObs = datos.iloc[:,3:].count().sum()
        data_set1 = f'      {MaxConcObs}      {MaxFluxObs}      {MaxFluxCells} # Data Set 1: MaxConcObs, MaxFluxObs, MaxFluxCells\n'
        escritura.write(data_set1)
        
        data_set2 = f'{outnam}      {inConcObs}      {inFluxObs}      {inSaveObs} # Data Set 2: OUTNAM, inConcObs, inFluxObs, inSaveObs\n'
        escritura.write(data_set2)
        
        nConcObs = datos.iloc[:,3:].count().sum()
        data_set3= f'      {nConcObs}    {Cscale:E}      {iOutCobs}       {iConcLOG}       {iConcINTP} # Data Set 3: nConcObs, CScale, iOutCobs, iConcLOG, iConcINTP\n'
        escritura.write(data_set3)
        #escritura.write('# Data Set 4: COBSNAM, Layer, Row, Column, iComp, TimeObs, Roff, Coff, weight, COBS\n')
        
        for fila in range(datos.shape[0]):

            COBS_nam = str(datos['ID'][fila])
            row, column, r_off, c_off = coord_params(datos['X'][fila], datos['Y'][fila],
                                                 x0, y0, dx, dy)
            n_obs = 1
            for columna in range(3,datos.shape[1]):
                if datos.notnull().iloc[fila, columna]:
                    yr = list(datos.columns)[columna]
                    timeObs = (yr - year_ini) * 3600 * 24 * 365
                
                    datos_geom = geometria[(geometria['COLUMN']==column) & (geometria['ROW']==row)]
                    dz = datos_geom.iloc[0,4]-datos_geom.iloc[0,5]
                
                    COBS = datos.iloc[fila, columna] * dx * dy * dz * rho  # mg/l *(1g/1000mg)(1000l/1m3) *(m3)* (adim) = g
                    COBSNAM = COBS_nam + '_' + str(n_obs)
                    n_obs += 1
                    data_set4 = f'{COBSNAM}     {layer}     {row}     {column}     {iComp}  {timeObs:E}   {r_off:E}   {c_off:E}   {weight:E}   {COBS:E}\n'#datos del pozo
                    escritura.write(data_set4) 
                
        if MaxFluxObs > 0:
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
     



def piezometria2ob_hob(outnam='output', rename=False, x0=455204.440, y0=2110063.17+64000, 
                       dx=2000, dy=2000, ml_obs=0, max_m=2, iu_hobsv=42, 
                       hob_dry=1.0E+30, tm_of_mult_hbs=1.0, 
                       na_val = ['ND','nd'], layer=2, year=1934, toffset = 3600 * 24 * 365):
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
        
        if rename:
            reg_rename = open(Path.cwd()/ 'ob_hob'/'salida'/'cambio_nombres.txt', 'w')
            reg_rename.write('Anterior\tActual\n')
        obs = datos.iloc[:,3:]
        numero_hobs = obs.count().sum()
        obs_bool = obs.isnull()
        archivo = outnam + '.ob_hob'
        encabezado = f'#Archivo ob_hob\t Created on {get_time()} by Dpto Rec. Nat, IGF UNAM \n'
        escritura = open(Path.cwd()/ 'ob_hob'/'salida'/archivo, 'w')
        escritura.write(encabezado)         
        
        dat_set1 = f'  {numero_hobs}     {ml_obs}    {max_m}    {iu_hobsv}   {hob_dry:E} # Data Set 1: NH MOBS MAXM IUHOBSV HOBDRY\n'
        escritura.write(dat_set1)
        
        dat_set2 = f"{  tm_of_mult_hbs:.6f}#Data Set 2: TOMULTH\n"
        escritura.write(dat_set2)       
        
        
        for fila in range(obs.shape[0]):
            obsname = str(datos['ID'][fila])
            if rename:
                anterior, obsname = corregir_nombre(obsname)
                reg_rename.write(f'{anterior}\t\t{obsname}\n')
                
                
            row, column, r_off, c_off = coord_params(datos['X'][fila], datos['Y'][fila], x0, y0, dx, dy)
            irefsp = -obs.loc[fila].count()
            itt = 2 #1 para usar observaciones de carga, 2 para usar cambios de carga como observaciones
            obs_num = 1
            for columna in range(obs.shape[1]):
                if not obs_bool.iloc[fila,columna]:
                    hobs = obs.iloc[fila,columna]
                    if obs_num == 1:
                        dat_set3 = f'{obsname}    {layer}    {row}    {column}   {irefsp}   {toffset:E}   {r_off:E}   {c_off:E}   {hobs:E} # Data Set 3: OBSNAM LAYER ROW COLUMN IREFSP TOFFSET ROFF COFF HOBS\n'
                        escritura.write(dat_set3)
                        dat_set5 = f'     {itt} # Data Set 5: ITT\n'
                        escritura.write(dat_set5)
                    yr = list(obs.columns)[columna]
                    obs_subname =  obsname + '_' + str(obs_num)
                    obs_num += 1
                    irefsp = int((yr) - year) + 1
                    dat_set6 = f'{obs_subname}    {irefsp}   {toffset:E}   {hobs:E} # Data Set 6: OBSNAM IREFSP TOFFSET HOBS\n'
                    escritura.write(dat_set6)
        if rename:
            reg_rename.close()
        escritura.close()
        print(f"Finalizado: archivo exportado en {Path.cwd()/ 'ob_hob'/'salida'/archivo}")
    else:
        print('Error al cargar el archivo de datos, verifique que se encuentre el archivo correctamente')
        print(f'En el directorio de entrada se encuentran {len(lista_piezometria)} archivos, podría estar abrierto el archivo')



def mk_pfile(outname='pesos', datos=None, sum_sd=None, sum_sd_dr=None, rename = False):
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
        if rename:
            n_nombres = open(Path.cwd()/ 'pesos'/'salida'/'cambio_nombres.txt', 'w')
            n_nombres.write('Anterior\tActual\n')        
        archivo = open(Path.cwd()/ 'pesos'/'salida'/nombre, 'w')
        for fila in range(datos.shape[0]):
            pozo = datos['ID'][fila]
            if rename:
                anterior, pozo = corregir_nombre(pozo)
                n_nombres.write(f'{anterior}\t\t{pozo}\n')
                
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
        if rename:
            n_nombres.close()
        archivo.close()
        print(f"Finalizado: archivo exportado en {Path.cwd()/ 'pesos'/'salida'/nombre}")
        
        
def pesos(nombre='pesos', calcular=False, na_val = ['ND','nd'], rename=False):
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
                    
            sd3.fillna(value=sd3['sd3'].mean(), inplace=True)
            sd3.set_index('ID', inplace = True)
            #print(sd3)

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
            #print(sd4)
            
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
            
            mk_pfile(nombre,hojas[0], sum_sd_carga,sum_sd_abat, rename=rename)
            
            #Exportar los datos a un excel para volver a crear el archivo de pesos
            writer = pd.ExcelWriter(Path.cwd()/ 'pesos'/'salida'/'datos_completo.xlsx', engine='xlsxwriter')
            hojas[0].fillna(value='ND').to_excel(writer, sheet_name='datos')
            sd1.to_excel(writer, sheet_name='sd1')
            sd2.to_excel(writer, sheet_name='sd2')
            sd3.to_excel(writer, sheet_name='sd3')
            sd4.to_excel(writer, sheet_name='sd4')
            sd5.to_excel(writer, sheet_name='sd5')
            sum_sd_carga.to_excel(writer, sheet_name='suma_sd')
            sum_sd_abat.to_excel(writer, sheet_name='suma_sd_drawdowns')
            writer.save()
            print(f"Archivo de Excel exportado en {Path.cwd()/'pesos'/'salida'/'datos_completo.xlsx'}")
        
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
            mk_pfile(nombre, datos, suma_sd, suma_sd_dr, rename=rename)

