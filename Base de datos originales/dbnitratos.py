# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:04:14 2020

@author: Jorge Antonio Matías López
"""
import numpy as np
import pandas as pd
from pathlib import Path


def ls(ruta = Path.cwd()):
    '''Función que sirve para obtener una lista con los nombres
       de los archivos contenidos en una carpeta 
    '''
    return [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]

nombres = ls()
conjunto_db = []

for nombre in nombres:
    df = pd.read_excel(Path.cwd() / nombre)
    conjunto_db.append(df)

salida = Path.cwd() / 'output' # ruta a la carpeta de salida para almacenar las gráficas

#obsformation = pd.read_csv(Path.cwd() / 'input' / 'obs_formation.csv') #lectura del csv obs_formation
    
#df = pd.read_fwf(Path.cwd().joinpath('input').joinpath('.hob_out'))