# -*- coding: utf-8 -*-
"""
Created on Thu Aug 1 2020

@author: Jorge Antonio Matías López
"""
from pathlib import Path

def contar_archivo(ruta = Path.cwd(), ext='xls'):
    files = [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]
    lista = []
    for archivo in files:
        if archivo.count('xls') != 0:
            lista.append(archivo)
    return lista

    
