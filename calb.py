#!/usr/bin/python
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import dates

import numpy as np
import pylab
import pandas as pd
import random
from datetime import datetime, timedelta
from collections import OrderedDict

from scipy.io.idl import readsav
from scipy.optimize import curve_fit
from scipy import stats
from scipy import signal
from IPython.core.display import Image
import io
import base64
from IPython.display import HTML

import sunpy
from astropy.io import fits
from sunpy.sun import constants as con
#from sunpy.net.helioviewer import HelioviewerClient
from sunpy.time import *
#from sunpy.net import vso
#from sunpy import lightcurve as lc
from sunpy.time import TimeRange
#from sunpy.net import hek

from attrdict import AttrDict
from astropy.constants import k_B,h

#import sys
#sys.path.append("/home/fer/Documents/myPYTHONstuff/modules/")
#from mapping import *
#from other import *

def hola():
    print('Hello World, this is my calibration module of SST raw data saved in IDL SAVE file type')
    return

def integrador(myarray,janela,tipo):
    forma = np.shape(myarray)
    valor = forma[0]//janela-1
    difer = forma[0]-valor*janela
    if len(forma) <= 1:
        aux = myarray[:-difer]
        if tipo =='time':
            matriz = np.reshape(aux,(valor,janela))
            salida = matriz[:,0][:]
        if tipo =='data':
            matriz = np.reshape(aux,(valor,janela))
            salida = np.sum(matriz,axis=1)/janela
    else:
        salida = np.zeros(shape=(forma[0]-difer,forma[1]))
        for k in range(forma[1]):
            aux = myarray[:-difer,k]
            matriz = np.reshape(aux,(valor,janela))
            previa = np.sum(matriz,axis=1)/janela
            salida[:,k] = previa
    return salida

def lee_sst(pathplusfile):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (bi or rs = biors).
    """
    biorsin = readsav(pathplusfile,python_dict=False,verbose=False)
    biorsout= AttrDict({'pm_daz':[],'off':[],'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],
                     'pos_time':[],'recnum':[],'opmode':[],'time':[],'pm_del':[],'gps_status':[],'adcval':[],
                     'y_off':[],'t     arget':[]})
    for items in biorsin:
        biorsout[items] = biorsin[items].copy()
    try:
        yyyy = np.int(pathplusfile[-11:-9])+2000
        mm = np.int(pathplusfile[-9:-7])
        dd = np.int(pathplusfile[-7:-5])
    except:
        yyyy = np.int(pathplusfile[-16:-14])+2000
        mm = np.int(pathplusfile[-14:-12])
        dd = np.int(pathplusfile[-12:-10])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    thatday = obsday + biorsout['time']/1.0e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biorsout

def obsolete_time(rs,yyyy,mm,dd):
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    thatday = obsday + rs['time']/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    ti = dates.date2num(dts)
    return ti,dts
'''
def lee2sst(fileplace1,fileplace2):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors) de dos archivos consecutivos y concatenados.
    """
    try:
	#pathplusfile1,pathplusfile2 = fileplace1,fileplace2
	lenstr = len(pathplusfile1)
	biors1 = readsav(pathplusfile1,python_dict=False,verbose=False)
	biors2 = readsav(pathplusfile2,python_dict=False,verbose=False)
	yyyy = np.int(pathplusfile1[lenstr-16:lenstr-14])+2000
	mm = np.int(pathplusfile1[lenstr-14:lenstr-12])
	dd = np.int(pathplusfile1[lenstr-12:lenstr-10])
        obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
        biors12= AttrDict({'pm_daz':[],'off':[],'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],'pos_time':[],
        'recnum':[],'opmode':[],'time':[],'pm_del':[],'gps_status':[],'adcval':[],'y_off':[],'target':[]})
        for items in biors12:
            a,b = biors1[items],biors2[items]
	    biors12[items] = np.concatenate((a,b))
        thatday = obsday + biors12['time']/1.e4
        dts = map(datetime.fromtimestamp,thatday)
        fds = dates.date2num(dts)
    except:
        print('no implementado aun')
    return fds,biors12
'''
def lee2oldsst(fileplace1,fileplace2):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors) de dos archivos consecutivos y concatenados.
    """

    pathplusfile1,pathplusfile2 = fileplace1,fileplace2
    lenstr = len(pathplusfile1)
    biors1 = readsav(pathplusfile1,python_dict=False,verbose=False)
    biors2 = readsav(pathplusfile2,python_dict=False,verbose=False)
    yyyy = np.int(pathplusfile1[lenstr-16:lenstr-14])+2000
    mm = np.int(pathplusfile1[lenstr-14:lenstr-12])
    dd = np.int(pathplusfile1[lenstr-12:lenstr-10])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    biors12= AttrDict({'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],'pos_time':[],
    'recnum':[],'opmode':[],'time':[],'gps_status':[],'adcval':[],'y_off':[],'target':[]})
    for items in biors1:
        a,b = biors1[items],biors2[items]
        biors12[items] = np.concatenate((a,b))
    thatday = obsday + biors12['time']/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biors12

def lee3sst(fileplace1,fileplace2,fileplace3):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors) de dos archivos consecutivos y concatenados.
    """

    pathplusfile1,pathplusfile2,pathplusfile3 = fileplace1,fileplace2,fileplace3
    lenstr = len(pathplusfile1)
    biors1 = readsav(pathplusfile1,python_dict=True,verbose=False)
    biors2 = readsav(pathplusfile2,python_dict=True,verbose=False)
    biors3 = readsav(pathplusfile3,python_dict=True,verbose=False)
    yyyy = np.int(pathplusfile1[lenstr-16:lenstr-14])+2000
    mm = np.int(pathplusfile1[lenstr-14:lenstr-12])
    dd = np.int(pathplusfile1[lenstr-12:lenstr-10])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    biors123= AttrDict({'pm_daz':[],'off':[],'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],'pos_time':[],
    'recnum':[],'opmode':[],'time':[],'pm_del':[],'gps_status':[],'adcval':[],'y_off':[],'target':[]})
    for items in biors1:
        a,b,c = biors1[items],biors2[items],biors3[items]
        biors123[items] = np.concatenate((a,b,c))
    thatday = obsday + biors123['time']/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biors123

def lee4sst(fileplace1,fileplace2,fileplace3,fileplace4):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors) de dos archivos consecutivos y concatenados.
    """
    pathplusfile1,pathplusfile2,pathplusfile3,pathplusfile4 = fileplace1,fileplace2,fileplace3,fileplace4
    lenstr = len(pathplusfile1)
    biors1 = readsav(pathplusfile1,python_dict=False,verbose=False)
    biors2 = readsav(pathplusfile2,python_dict=False,verbose=False)
    biors3 = readsav(pathplusfile3,python_dict=False,verbose=False)
    biors4 = readsav(pathplusfile4,python_dict=False,verbose=False)
    yyyy = np.int(pathplusfile1[lenstr-16:lenstr-14])+2000
    mm = np.int(pathplusfile1[lenstr-14:lenstr-12])
    dd = np.int(pathplusfile1[lenstr-12:lenstr-10])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    biors1234= AttrDict({'pm_daz':[],'off':[],'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],'pos_time':[],
    'recnum':[],'opmode':[],'time':[],'pm_del':[],'gps_status':[],'adcval':[],'y_off':[],'target':[]})
    for items in biors1:
        a,b,c,d = biors1[items],biors2[items],biors3[items],biors4[items]
        biors1234[items] = np.concatenate((a,b,c,d))
    thatday = obsday + biors1234['time']/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biors1234

def lee5sst(fileplace1,fileplace2,fileplace3,fileplace4,fileplace5):
    """
    Lee archivos '/path/to/file/bi1yymmdd.save' o
    '/path/to/file/rs1yymmdd.hh00.save' y devuelve un
    vector-tiempo (fds) y la estructura completa de da-
    tos (biors) de dos archivos consecutivos y concatenados.
    """

    pathplusfile1,pathplusfile2,pathplusfile3,pathplusfile4,pathplusfile5 = fileplace1,fileplace2,fileplace3,fileplace4,fileplace5
    lenstr = len(pathplusfile1)
    biors1 = readsav(pathplusfile1,python_dict=False,verbose=False)
    biors2 = readsav(pathplusfile2,python_dict=False,verbose=False)
    biors3 = readsav(pathplusfile3,python_dict=False,verbose=False)
    biors4 = readsav(pathplusfile4,python_dict=False,verbose=False)
    biors5 = readsav(pathplusfile5,python_dict=False,verbose=False)
    yyyy = np.int(pathplusfile1[lenstr-16:lenstr-14])+2000
    mm = np.int(pathplusfile1[lenstr-14:lenstr-12])
    dd = np.int(pathplusfile1[lenstr-12:lenstr-10])
    obsday = int(datetime(yyyy,mm,dd,0,0).strftime('%s'))
    biors12345= AttrDict({'pm_daz':[],'off':[],'azierr':[],'eleerr':[],'x_off':[],'elepos':[],'azipos':[],'pos_time':[],
    'recnum':[],'opmode':[],'time':[],'pm_del':[],'gps_status':[],'adcval':[],'y_off':[],'target':[]})
    for items in biors1:
        a,b,c,d,e = biors1[items],biors2[items],biors3[items],biors4[items],biors5[items]
        biors12345[items] = np.concatenate((a,b,c,d,e))
    thatday = obsday + biors12345['time']/1.e4
    dts = map(datetime.fromtimestamp,thatday)
    fds = dates.date2num(dts)
    return fds,biors12345

def modos(biors):
    """
    Devuelve en pantalla la lista de modos existentes
    en el archivo 'biors' y los guarda en el diccionario
    'opmodict'.
    """
    modos=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,40,50,55,99]
    opmodict={'0':0,'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0,'9':0,'10':0,'11':0,\
          '12':0,'13':0,'14':0,'15':0,'16':0,'17':0,'18':0,'19':0,'20':0,'21':0,'22':0,\
          '23':0,'40':0,'50':0,'55':0,'99':0}
    for i in modos:
        xx=np.where(biors['opmode'] == i)
        tama = np.size(xx)
        if (tama >  0):
            opmodict[str(i)] = tama
            print( i,tama)
    cc = np.where(biors['target']/32==1)
    hh = np.where(biors['target']/32==2)
    if (np.size(cc) > 0):
        print( 'cold source')
    if (np.size(hh) > 0):
        print( 'hot source')
    return opmodict

def taxobi(fdsh,bh):
    """
    Separa los datos adcval y los correspondientes
    vectores-tiempo segun los modos mas comunes del
    'opmode'. Guarda los datos segregados en el dic-
    cionario 'taxdict'. Funciona solo con archivos de
    tipo rs1yymmdd.hh00.save
    """
    #####################################################
    colfds = fdsh[np.where(bh['target']/32 == 1)]
    coladc = bh['adc'][np.where(bh['target']/32 == 1)]
    hotfds = fdsh[np.where(bh['target']/32 == 2)]
    hotadc = bh['adc'][np.where(bh['target']/32 == 2)]
    #####################################################
    trkadc = bh['adc'][np.where(bh['opmode'] == 0)]
    trkfds = fdsh[np.where(bh['opmode'] == 0)]
    rdmadc = bh['adc'][np.where(bh['opmode'] == 1)]
    rdmfds = fdsh[np.where(bh['opmode'] == 1)]
    mapadc = bh['adc'][np.where(bh['opmode'] == 2)]
    mapfds = fdsh[np.where(bh['opmode'] == 2)]
    rmpadc = bh['adc'][np.where(bh['opmode'] == 3)]
    rmpfds = fdsh[np.where(bh['opmode'] == 3)]
    mpiadc = bh['adc'][np.where(bh['opmode'] == 4)]
    mpifds = fdsh[np.where(bh['opmode'] == 4)]
    scnadc = bh['adc'][np.where(bh['opmode'] == 5)]
    scnfds = fdsh[np.where(bh['opmode'] == 5)]
    sciadc = bh['adc'][np.where(bh['opmode'] == 9)]
    scifds = fdsh[np.where(bh['opmode'] == 9)]
    tauadc = bh['adc'][np.where(bh['opmode'] == 10)]
    taufds = fdsh[np.where(bh['opmode'] == 10)]
    stladc = bh['adc'][np.where(bh['opmode'] == 50)]
    stlfds = fdsh[np.where(bh['opmode'] == 50)]
    unkadc = bh['adc'][np.where(bh['opmode'] == 99)]
    unkfds = fdsh[np.where(bh['opmode'] == 99)]
    taxdict={\
             'coltime':colfds,'colval':coladc\
             ,'hottime':hotfds,'hotval':hotadc\
             ,'trktime':trkfds,'trkval':trkadc\
             ,'rdmtime':rdmfds,'rdmval':rdmadc\
             ,'maptime':mapfds,'mapval':mapadc\
             ,'rmptime':rmpfds,'rmpval':rmpadc\
             ,'mpitime':mpifds,'mpival':mpiadc\
             ,'scntime':scnfds,'scnval':scnadc\
             ,'scitime':scifds,'scival':sciadc\
             ,'tautime':taufds,'tauval':tauadc\
             ,'stltime':stlfds,'stlval':stladc\
             ,'unktime':unkfds,'unkval':unkadc\
             }
    print(np.size(taxdict['colval']),np.size(taxdict['hotval']))
    return taxdict

def taxonomia(fdsh,bh):
    """
    Separa los datos adcval y los correspondientes
    vectores-tiempo segun los modos mas comunes del
    'opmode'. Guarda los datos segregados en el dic-
    cionario 'taxdict'. Funciona solo con archivos de
    tipo rs1yymmdd.hh00.save
    """
    #####################################################
    colfds = fdsh[np.where(bh['target']/32 == 1)]
    coladc = bh['adcval'][np.where(bh['target']/32 == 1)]
    hotfds = fdsh[np.where(bh['target']/32 == 2)]
    hotadc = bh['adcval'][np.where(bh['target']/32 == 2)]
    #####################################################
    trkadc = bh['adcval'][np.where(bh['opmode'] == 0)]
    trkfds = fdsh[np.where(bh['opmode'] == 0)]
    rdmadc = bh['adcval'][np.where(bh['opmode'] == 1)]
    rdmfds = fdsh[np.where(bh['opmode'] == 1)]
    mapadc = bh['adcval'][np.where(bh['opmode'] == 2)]
    mapfds = fdsh[np.where(bh['opmode'] == 2)]
    rmpadc = bh['adcval'][np.where(bh['opmode'] == 3)]
    rmpfds = fdsh[np.where(bh['opmode'] == 3)]
    mpiadc = bh['adcval'][np.where(bh['opmode'] == 4)]
    mpifds = fdsh[np.where(bh['opmode'] == 4)]
    scnadc = bh['adcval'][np.where(bh['opmode'] == 5)]
    scnfds = fdsh[np.where(bh['opmode'] == 5)]
    sciadc = bh['adcval'][np.where(bh['opmode'] == 9)]
    scifds = fdsh[np.where(bh['opmode'] == 9)]
    tauadc = bh['adcval'][np.where(bh['opmode'] == 10)]
    taufds = fdsh[np.where(bh['opmode'] == 10)]
    stladc = bh['adcval'][np.where(bh['opmode'] == 50)]
    stlfds = fdsh[np.where(bh['opmode'] == 50)]
    unkadc = bh['adcval'][np.where(bh['opmode'] == 99)]
    unkfds = fdsh[np.where(bh['opmode'] == 99)]
    taxdict={\
             'coltime':colfds,'colval':coladc\
             ,'hottime':hotfds,'hotval':hotadc\
             ,'trktime':trkfds,'trkval':trkadc\
             ,'rdmtime':rdmfds,'rdmval':rdmadc\
             ,'maptime':mapfds,'mapval':mapadc\
             ,'rmptime':rmpfds,'rmpval':rmpadc\
             ,'mpitime':mpifds,'mpival':mpiadc\
             ,'scntime':scnfds,'scnval':scnadc\
             ,'scitime':scifds,'scival':sciadc\
             ,'tautime':taufds,'tauval':tauadc\
             ,'stltime':stlfds,'stlval':stladc\
             ,'unktime':unkfds,'unkval':unkadc\
             }
    print(np.size(taxdict['colval']),np.size(taxdict['hotval']))
    return taxdict

def colecta(taxdict):
    """
    Escoge percentiles varios  de cada modo encontrado en
    'taxonomia' y compone un diccionario de
    estadisticas 'codict'.
    """
    codict = {}
    if len(taxdict['colval']) > 0:
        coldadc = taxdict['colval']
        statcol = [np.percentile(coldadc,12,axis=0),np.percentile(coldadc,25,axis=0),
                np.percentile(coldadc,37,axis=0),np.percentile(coldadc,50,axis=0),
                np.percentile(coldadc,62,axis=0),np.percentile(coldadc,65,axis=0)]
        extra = {'colstat':statcol}
        codict.update(extra)
    if len(taxdict['hotval']) > 0:
        hotadc = taxdict['hotval']
        stathot = [np.percentile(hotadc,12,axis=0),np.percentile(hotadc,25,axis=0),
                np.percentile(hotadc,37,axis=0),np.percentile(hotadc,50,axis=0),
                np.percentile(hotadc,62,axis=0),np.percentile(hotadc,65,axis=0)]
        extra = {'hotstat':stathot}
        codict.update(extra)

    if len(taxdict['trkval']) > 0:
        trkadc = taxdict['trkval']
        stattrk = [np.percentile(trkadc,12,axis=0),np.percentile(trkadc,23,axis=0),
                np.percentile(trkadc,30,axis=0),np.percentile(trkadc,50,axis=0),
                np.percentile(trkadc,60,axis=0),np.percentile(trkadc,65,axis=0)]
        extra = {'trkstat':stattrk}
        codict.update(extra)
    if len(taxdict['mapval']) > 0:
        mapadc = taxdict['mapval']
        statmap = [np.percentile(mapadc,12,axis=0),np.percentile(mapadc,23,axis=0),
                np.percentile(mapadc,30,axis=0),np.percentile(mapadc,50,axis=0),
                np.percentile(mapadc,60,axis=0),np.percentile(mapadc,65,axis=0)]
        extra = {'mapstat':statmap}
        codict.update(extra)
    if len(taxdict['mpival']) > 0:
        mpiadc = taxdict['mpival']
        statmpi = [np.percentile(mpiadc,12,axis=0),np.percentile(mpiadc,25,axis=0),
                np.percentile(mpiadc,37,axis=0),np.percentile(mpiadc,50,axis=0),
                np.percentile(mpiadc,62,axis=0),np.percentile(mpiadc,65,axis=0)]
        extra = {'mpistat':statmpi}
        codict.update(extra)
    if len(taxdict['scnval']) > 0:
        scnadc = taxdict['scnval']
        statscn = [np.percentile(scnadc,12,axis=0),np.percentile(scnadc,23,axis=0),
                np.percentile(scnadc,30,axis=0),np.percentile(scnadc,50,axis=0),
                np.percentile(scnadc,60,axis=0),np.percentile(scnadc,65,axis=0)]
        extra = {'scnstat':statscn}
        codict.update(extra)
    if len(taxdict['scival']) > 0:
        sciadc = taxdict['scival']
        statsci = [np.percentile(sciadc,12,axis=0),np.percentile(sciadc,23,axis=0),
                np.percentile(sciadc,30,axis=0),np.percentile(sciadc,50,axis=0),
                np.percentile(sciadc,60,axis=0),np.percentile(sciadc,65,axis=0)]
        extra = {'scistat':statsci}
        codict.update(extra)
    if len(taxdict['tauval']) > 0:
        tauadc = taxdict['tauval']
        stattau = [np.percentile(tauadc,12,axis=0),np.percentile(tauadc,25,axis=0),
                np.percentile(tauadc,37,axis=0),np.percentile(tauadc,50,axis=0),
                np.percentile(tauadc,62,axis=0),np.percentile(tauadc,65,axis=0)]
        extra = {'taustat':stattau}
        codict.update(extra)
    if len(taxdict['stlval']) > 0:
        stladc = taxdict['stlval']
        statstl = [np.percentile(stladc,12,axis=0),np.percentile(stladc,25,axis=0),
                np.percentile(stladc,37,axis=0),np.percentile(stladc,50,axis=0),
                np.percentile(stladc,62,axis=0),np.percentile(stladc,65,axis=0)]
        extra = {'stlstat':statstl}
        codict.update(extra)
    """
    if len(taxdict['unkval']) > 0:
        unkadc = taxdict['unkval']
        statunk = [\
        np.percentile(unkadc,25,axis=0),np.percentile(unkadc,50,axis=0),np.percentile(unkadc,75,axis=0)]
        extra = {'unkstat':statunk}
        codict.update(extra)
    """
    return codict

def preparacal(codict):
    """
    Se trata de colocar las estadisticas de obtenidas de
    'colecta' y almacenados en el dicccionario 'codict'
    para hacer la pseudo-calibracion que iguala canales.
    """
    prdict = {}
    stat1 = [0.0]
    stat2 = [0.0]
    stat3 = [0.0]
    stat4 = [0.0]
    stat5 = [0.0]
    stat6 = [0.0]
    for key in codict.keys():
        var1 = [codict[key][0][0],codict[key][1][0],codict[key][2][0],codict[key][3][0],codict[key][4][0],codict[key][5][0]]
        stat1.extend(var1)
        var2 = [codict[key][0][1],codict[key][1][1],codict[key][2][1],codict[key][3][1],codict[key][4][1],codict[key][5][1]]
        stat2.extend(var2)
        var3 = [codict[key][0][2],codict[key][1][2],codict[key][2][2],codict[key][3][2],codict[key][4][2],codict[key][5][2]]
        stat3.extend(var3)
        var4 = [codict[key][0][3],codict[key][1][3],codict[key][2][3],codict[key][3][3],codict[key][4][3],codict[key][5][3]]
        stat4.extend(var4)
        var5 = [codict[key][0][4],codict[key][1][4],codict[key][2][4],codict[key][3][4],codict[key][4][4],codict[key][5][4]]
        stat5.extend(var5)
        var6 = [codict[key][0][5],codict[key][1][5],codict[key][2][5],codict[key][3][5],codict[key][4][5],codict[key][5][5]]
        stat6.extend(var6)
    del stat1[0],stat2[0],stat3[0],stat4[0],stat5[0],stat6[0]
    prdict = {'stat01':stat1,'stat02':stat2,'stat03':stat3,'stat04':stat4,'stat05':stat5,'stat06':stat6}
    return prdict

def pseudocalbi(prdict,rsHH):
    '''
    La tan esperada y poco entendida pseudocalibracion.
    Tomo los valores de rsHH y los alineo segun la com-
    paracion con la estadistica derivada de prdict
    (diccionario obtenido con preparacal)
    '''
    #CANALES 212
    degree = 2
    #eixo X
    rsHH['adc'][:,0] -= np.nanmin(rsHH['adc'][:,0])
    xi = prdict['stat01']
    A = np.vander(xi,degree)
    #eixos Y
    y = prdict['stat02']
    (coeffs2, residuals2, rank2, sing_vals2) = np.linalg.lstsq(A, y, rcond=None)
    #a,b,c = coeffs2[0],coeffs2[1],coeffs2[2]
    print(coeffs2)
    f2 = np.poly1d(coeffs2)
    rsHH['adc'][:,1] = (rsHH['adc'][:,1]-coeffs2[1])/coeffs2[0]
    rsHH['adc'][:,1] -=np.nanmin(rsHH['adc'][:,1])
    #integrato12.adcval[:,1] = np.sqrt((integrato12.adcval[:,1]-c)/a + b**2/(4*a))-b/(2*a)
    y = prdict['stat03']
    (coeffs3, residuals3, rank3, sing_vals3) = np.linalg.lstsq(A, y, rcond=None)
    #a,b,c = coeffs3[0],coeffs3[1],coeffs3[2]
    print(coeffs3)
    f3 = np.poly1d(coeffs3)
    rsHH['adc'][:,2] = (rsHH['adc'][:,2]-coeffs3[1])/coeffs3[0]
    rsHH['adc'][:,2] -=np.nanmin(rsHH['adc'][:,2])
    #integrato12.adcval[:,2] = np.sqrt((integrato12.adcval[:,2]-c)/a + b**2/(4*a))-b/(2*a)
    y = prdict['stat04']
    (coeffs4, residuals4, rank4, sing_vals4) = np.linalg.lstsq(A, y, rcond=None)
    #a,b,c = coeffs4[0],coeffs4[1],coeffs4[2]
    print(coeffs4)
    f4 = np.poly1d(coeffs4)
    rsHH['adc'][:,3] = (rsHH['adc'][:,3]-coeffs4[1])/coeffs4[0]
    rsHH['adc'][:,3] -=np.nanmin(rsHH['adc'][:,3])
    #integrato12.adcval[:,3] = np.sqrt((integrato12.adcval[:,3]-c)/a + b**2/(4*a))-b/(2*a
    
    #CANALES 405
    degree = 2
    #eixo X
    rsHH['adc'][:,5] -= np.nanmin(rsHH['adc'][:,0]) 
    xi = prdict['stat06']
    A = np.vander(xi,degree)
    #eixos Y
    y = prdict['stat05']
    (coeffs5, residuals5, rank5, sing_vals5) = np.linalg.lstsq(A, y, rcond=None)
    print(coeffs5)
    f5 = np.poly1d(coeffs5)
    rsHH['adc'][:,4] = (rsHH['adc'][:,4]-coeffs5[1])/coeffs5[0]
    rsHH['adc'][:,4] -=np.nanmin(rsHH['adc'][:,4])
    return rsHH

def pseudocal(prdict,rsHH):
    '''
    La tan esperada y poco entendida pseudocalibracion.
    Tomo los valores de rsHH y los alineo segun la com-
    paracion con la estadistica derivada de prdict
    (diccionario obtenido con preparacal)
    '''
    #CANALES 212
    degree = 2

    #eixo X

    rsHH['adcval'][:,0] -= np.nanmin(rsHH['adcval'][:,0])
    xi = prdict['stat01']
    A = np.vander(xi,degree)
    
    #eixos Y
    
    y = prdict['stat02']
    (coeffs2, residuals2, rank2, sing_vals2) = np.linalg.lstsq(A, y, rcond=None)
    print(coeffs2)
    f2 = np.poly1d(coeffs2)
    aux2  = rsHH['adcval'][:,1].astype(float)
    aux2 = (aux2-coeffs2[1])/coeffs2[0]
    aux2 -=np.nanmin(aux2)
    rsHH['adcval'][:,1] = aux2
    
    y = prdict['stat03']
    (coeffs3, residuals3, rank3, sing_vals3) = np.linalg.lstsq(A, y, rcond=None)
    print(coeffs3)
    f3 = np.poly1d(coeffs3)
    aux3  = rsHH['adcval'][:,2].astype(float)
    aux3 = (aux3-coeffs3[1])/coeffs3[0]
    aux3 -=np.nanmin(aux3)
    rsHH['adcval'][:,2] = aux3

    y = prdict['stat04']
    (coeffs4, residuals4, rank4, sing_vals4) = np.linalg.lstsq(A, y, rcond=None)
    print(coeffs4)
    f4 = np.poly1d(coeffs4)
    aux4  = rsHH['adcval'][:,3].astype(float)
    aux4 = (aux4-coeffs4[1])/coeffs4[0]
    aux4 -=np.nanmin(aux4)
    rsHH['adcval'][:,3] = aux4
    
    #CANALES 405
    degree = 2
    
    #eixo X
    
    rsHH['adcval'][:,5] -= np.nanmin(rsHH['adcval'][:,0]) 
    xi = prdict['stat06']
    A = np.vander(xi,degree)
    
    #eixos Y
    
    y = prdict['stat05']
    (coeffs5, residuals5, rank5, sing_vals5) = np.linalg.lstsq(A, y, rcond=None)
    print(coeffs5)
    f5 = np.poly1d(coeffs5)
    aux5  = rsHH['adcval'][:,4].astype(float)
    aux5 = (aux5-coeffs5[1])/coeffs5[0]
    aux5 -=np.nanmin(aux5)
    rsHH['adcval'][:,4] = aux5
    return rsHH

def tempcalbi(iicoHH,rsHH,t1=15.,t2_212=150.,t2_405=150.,ice=[1,1,1,1,1,1]):
    '''
    Calibracion en temperatura, toma los valores de la
    temperatura ambiente y la fuente caliente de los
    datos rs, necesita tambien los respectivos valores
    de ADC, colectados en iicoHH para finalmente corre-
    gir los valores rsHH
    '''
    cool_list = np.array(iicoHH['colstat'])
    #print cool_list
    hot_list = np.array(iicoHH['hotstat'])
    #print hot_list
    min
    #####212#####
    for i in range(0,4):
        #print i
        ADC2_ADC1 = hot_list[1][i]-cool_list[1][i]
        #print np.mean(hot_list[:][i]),np.mean(cool_list[:][i])
        t1xADC2 = t1*hot_list[1][i]
        t2xADC1 = ice[i]*t2_212*cool_list[1][i]
        m212 = (ice[i]*t2_212-t1)/ADC2_ADC1
        b212 = (t1xADC2-t2xADC1)/ADC2_ADC1
        print(1./m212,b212)
        tempval  = 1.0*rsHH['adc'][:,i]*m212+b212
        minival = np.nanmin(tempval)
        rsHH['adc'][:,i] = tempval#-minival#np.where(tempval>0,tempval,0)
    #####405#####
    for i in range(4,6):
        #print i
        ADC2_ADC1 = hot_list[1][i]-cool_list[1][i]
        t1xADC2 = t1*hot_list[1][i]
        t2xADC1 = ice[i]*t2_405*cool_list[1][i]
        m405 = (ice[i]*t2_405-t1)/ADC2_ADC1
        b405 = (t1xADC2-t2xADC1)/ADC2_ADC1
        #for i  in range(4,6):
        print(1./m405,b405)
        otempval  = 1.0*rsHH['adc'][:,i]*m405+b405
        rsHH['adc'][:,i] = otempval####np.where(otempval>0,otempval,0)#rsHH['adcval'][:,i]*m405[i]+b405[i]
    return rsHH

def tempcal(iicoHH,rsHH,t1=15.,t2_212=150.,t2_405=150.,ice=[1,1,1,1,1,1],tbase=60):
    '''
    Calibracion en temperatura, toma los valores de la
    temperatura ambiente y la fuente caliente de los
    datos rs, necesita tambien los respectivos valores
    de ADC, colectados en iicoHH para finalmente corre-
    gir los valores rsHH
    '''
    ###############################################
    ##scantauori=rsHH['adcval'][rsHH['opmode']==10,:]
    ##scanintori=rsHH['adcval'][rsHH['opmode']==9,:]
    ############################################### 

    cool_list = np.array(iicoHH['colstat'])
    #print cool_list
    hot_list = np.array(iicoHH['hotstat'])
    #print hot_list
    min
    #####212#####
    for i in range(0,4):
        ADC2_ADC1 = hot_list[1][i]-cool_list[1][i]
        #print np.mean(hot_list[:][i]),np.mean(cool_list[:][i])
        t1xADC2 = t1*hot_list[1][i]
        t2xADC1 = ice[i]*t2_212*cool_list[1][i]
        m212 = (ice[i]*t2_212-t1)/ADC2_ADC1
        b212 = (t1xADC2-t2xADC1)/ADC2_ADC1
        print(1./m212,b212)
        tempval  = 1.0*rsHH['adcval'][:,i]*m212+b212
        #minval = np.nanmin(tempval)
        rsHH['adcval'][:,i] = tempval###-minval+tbase####np.where(tempval>0,tempval,0)
    #####405#####
    for i in range(4,6):
        ADC2_ADC1 = hot_list[1][i]-cool_list[1][i]
        t1xADC2 = t1*hot_list[1][i]
        t2xADC1 = ice[i]*t2_405*cool_list[1][i]
        m405 = (ice[i]*t2_405-t1)/ADC2_ADC1
        b405 = (t1xADC2-t2xADC1)/ADC2_ADC1
        #for i  in range(4,6):
        print(1./m405,b405)
        otempval  = 1.0*rsHH['adcval'][:,i]*m405+b405
        #ominval = np.nanmin(otempval)
        rsHH['adcval'][:,i] = otempval###-ominval+tbase####np.where(otempval>0,otempval,0)#rsHH['adcval'][:,i]*m405[i]+b405[i]
    return rsHH

def tempcal2(iicoHH,rsHH,t1=15.,t2_212=150.,t2_405=150.,ice=[0.999,0.999,0.999,0.999,0.999,0.999]):
    '''
    Calibracion en temperatura, toma los valores de la
    temperatura ambiente y la fuente caliente de los
    datos rs, necesita tambien los respectivos valores
    de ADC, colectados en iicoHH para finalmente corre-
    gir los valores rsHH
    '''
    zerado = np.nanmin([np.nanmin(rsHH['adcval'][rsHH['opmode']==5]),\
    #np.nanmin(rsHH['adcval'][rsHH['opmode']==5]),\
    #np.nanmin(rsHH['adcval'][rsHH['opmode']==9]),\
    #np.nanmin(rsHH['adcval'][rsHH['opmode']==50]),\
    #np.nanmin(rsHH['adcval'][rsHH['opmode']==99]),\
    ])

    rsHH['adcval'][rsHH['opmode']==99] = zerado
    rsHH['adcval'][rsHH['opmode']==50] = zerado
    rsHH['adcval'][rsHH['opmode']==10] = zerado
    rsHH['adcval'][rsHH['opmode']==9] = zerado
    myKnum = np.ones_like(rsHH['adcval'],dtype='float')-1.0
    myKden = np.ones_like(rsHH['adcval'],dtype='float')-1.0
    myK = np.ones_like(rsHH['adcval'],dtype='float')-1.0
    # = np.ones_like(rsHH)-1
    #cold_mean = np.array(iicoHH['colstat']).mean(axis=0)
    #hot_mean = np.array(iicoHH['hotstat']).mean(axis=0)
    #####212#####
    for i in range(0,4):
        cold_mean = np.array(iicoHH['colstat']).mean(axis=0)[i]
        hot_mean = np.array(iicoHH['hotstat']).mean(axis=0)[i]
        ##scn_mini = np.array(iicoHH['scnstat']).mean(axis=0)[i]
        #print scn_mini
        myKnum[:,i] = rsHH['adcval'][:,i] - cold_mean
        myKnum[:,i]-= np.min(myKnum[:,i])
        myKnum[:,i].clip(0.0)
        #if (np.nanmin(myKnum[:,i]) < 0):
        #    myKnum[:,i] -= np.min(myKnum[:,i])
        myKden[:,i] = hot_mean-cold_mean
        myK[:,i] = myKnum[:,i]/myKden[:,i]
        tempval  = myK[:,i]*(t2_212*ice[i]-t1)+t1
        rsHH['adcval'][:,i] = tempval - np.nanmin(tempval)#np.where(tempval>0,tempval,0)
        #rsHH['adcval'][:,i].clip(scn_mini,1.e6)
    #####405#####
    for i in range(4,6):
        cold_mean = np.array(iicoHH['colstat']).mean(axis=0)[i]
        hot_mean = np.array(iicoHH['hotstat']).mean(axis=0)[i]
        ##oscn_mini = np.array(iicoHH['scnstat']).mean(axis=0)[i]
        #print oscn_mini
        myKnum[:,i] = rsHH['adcval'][:,i] - cold_mean
        myKnum[:,i]-= np.min(myKnum[:,i])
        myKnum[:,i].clip(0.0)
        #if (np.nanmin(myKnum[:,i]) < 0):
        #    myKnum[:,i] -= np.min(myKnum[:,i])
        myKden[:,i] = hot_mean-cold_mean
        myK[:,i] = myKnum[:,i]/myKden[:,i]
        otempval  = myK[:,i]*(t2_405*ice[i]-t1)+t1
        rsHH['adcval'][:,i] = otempval - np.nanmin(otempval)#np.where(oempval>0,otempval,0)
        #rsHH['adcval'][:,i].clip(0.,1.e6)
    return rsHH,myKnum,myKden

def proxy_tau(iiitacc,rscc,eta_b212=0.17,eta_b405=0.09):
    '''
    Variables de entrada:
    iiitacc,rscc,eita212=0.32,eita405=0.16
    '''
    std_temp212,std_temp405 = np.std(iiitacc['scnval'][:,0:4],axis=0),np.std(iiitacc['scnval'][:,4:6],axis=0)
    mean_temp212,mean_temp405 = np.mean(iiitacc['scnval'][:,0:4],axis=0),np.mean(iiitacc['scnval'][:,4:6],axis=0)
    min_temp212,min_temp405 = np.min(iiitacc['scnval'][:,0:4],axis=0),np.min(iiitacc['scnval'][:,4:6],axis=0)
    max_temp212,max_temp405 = np.max(iiitacc['scnval'][:,0:4],axis=0),np.max(iiitacc['scnval'][:,4:6],axis=0)
    mean_angle, std_angle =  np.pi*np.mean(rscc['elepos']/1000.)/180.,np.pi*np.std(rscc['elepos']/1000.)/180.
    tauprox212 = np.sin(mean_angle)*np.log((5900.*eta_b212 - min_temp212)/(max_temp212-min_temp212))
    tauprox405 = np.sin(mean_angle)*np.log((5100.*eta_b405 - min_temp405)/(max_temp405-min_temp405))
    #print tauprox212,'\n',tauprox405
    return tauprox212,tauprox405

def opacalcbi(iiitaHH,rsHH,tbkg=30.,tau212=0.33,tsky212=280.,tau405=1.67,tsky405=265.):
    '''
    Calculo de opacidades a partir de valores calibrados,
    luego de la tercera taxonomia, siguiendo el modelo:
    tbkg*np.exp(-tau/x)+tsky*(1-np.exp(-tau/x))
    '''
    ###################theskytaufunction##################
    def skytaufunc(x,tbkg,tsky,tau):
        return tbkg*np.exp(-tau/x)+tsky*(1-np.exp(-tau/x))
    ######################################################
    rsHH['elepos'][rsHH['elepos']==0]=1.
    dtor = np.pi/180.
    x = np.sin(dtor*rsHH['elepos'][np.where(rsHH['opmode'] == 10)]/1.) # 1/airmass
    result212 = np.zeros((4,3))
    for i in range(0,4):
        y212 = iiitaHH['tauval'][:,i]
        inde = np.linspace(0,len(y212)-1,len(y212))
        z = np.polyfit([inde[0],inde[-1]],[y212[0],y212[-1]], 1)
        p = np.poly1d(z)
        y212 = y212-p(inde)
        try:
            popt212, pcov212 = curve_fit(skytaufunc, x, y212, [tbkg,tsky212,tau212],maxfev=600)
            result212[i,0],result212[i,1],result212[i,2]= popt212[2],popt212[1],popt212[0]
            print('ch',i+1,result212[i,0],result212[i,1],result212[i,2])
        except RuntimeError:
            print('ch',i+1,'Erro  no ajuste')
    result405 = np.zeros((2,3))
    for i in range(4,6):
        y405 = iiitaHH['tauval'][:,i]
        inde = np.linspace(0,len(y405)-1,len(y405))
        z = np.polyfit([inde[0],inde[-1]],[y405[0],y405[-1]],1)
        p = np.poly1d(z)
        y405 = y405-p(inde)
        try:
            popt405,pcov405 =  curve_fit(skytaufunc, x, y405, [tbkg,tsky405,tau405],maxfev=600)
            result405[i-4,0],result405[i-4,1],result405[i-4,2]= popt405[2],popt405[1],popt405[0]
            print('ch',i+1,result405[i-4,0],result405[i-4,1],result405[i-4,2])
        except RuntimeError:
            print('ch',i+1,'Erro no ajuste')
    return result212,result405

def opacalc(iiitaHH,rsHH,tbkg=30.,tau212=0.33,tsky212=280.,tau405=1.67,tsky405=265.):
    '''
    Calculo de opacidades a partir de valores calibrados,
    luego de la tercera taxonomia, siguiendo el modelo:
    tbkg*np.exp(-tau/x)+tsky*(1-np.exp(-tau/x))
    '''
    ###################theskytaufunction##################
    def skytaufunc(x,tbkg,tsky,tau):
        return tbkg*np.exp(-tau/x)+tsky*(1-np.exp(-tau/x))
    ######################################################
    rsHH['elepos'][rsHH['elepos']==0]=1000.
    dtor = np.pi/180.
    x = np.sin(dtor*rsHH['elepos'][np.where(rsHH['opmode'] == 10)]/1000.) # 1/airmass
    result212 = np.zeros((4,3))
    for i in range(0,4):
        y212 = iiitaHH['tauval'][:,i]
        inde = np.linspace(0,len(y212)-1,len(y212))
        z = np.polyfit([inde[0],inde[-1]],[y212[0],y212[-1]], 1)
        p = np.poly1d(z)
        y212 = y212-p(inde)
        try:
            popt212, pcov212 = curve_fit(skytaufunc, x, y212, [tbkg,tsky212,tau212],maxfev=600)
            result212[i,0],result212[i,1],result212[i,2]= popt212[2],popt212[1],popt212[0]
            print('ch',i+1,result212[i,0],result212[i,1],result212[i,2])
        except RuntimeError:
            print('ch',i+1,'Erro  no ajuste')
    result405 = np.zeros((2,3))
    for i in range(4,6):
        y405 = iiitaHH['tauval'][:,i]
        inde = np.linspace(0,len(y405)-1,len(y405))
        z = np.polyfit([inde[0],inde[-1]],[y405[0],y405[-1]],1)
        p = np.poly1d(z)
        y405 = y405-p(inde)
        try:
            popt405,pcov405 =  curve_fit(skytaufunc, x, y405, [tbkg,tsky405,tau405],maxfev=600)
            result405[i-4,0],result405[i-4,1],result405[i-4,2]= popt405[2],popt405[1],popt405[0]
            print('ch',i+1,result405[i-4,0],result405[i-4,1],result405[i-4,2])
        except RuntimeError:
            print('ch',i+1,'Erro no ajuste')
    return result212,result405

def opadescbi(rsHH,opa212=0.33,ts212=280.,opa405=1.67,ts405=265.0):
    '''
    Descuento de opacidad atmosferica, los tau son ingresados
    del site http://190.3.114.98:20081/sst/sst_plot.php/
    Requiere tambien la tercera taxonomia: iiitaHH (despues
    de la calibracion en temperatura) para modificar los datos
    rsHH, de lo cuales tambien toma las elevaciones. Tiene que
    existir calculo de opacidad en ese archivo rs.
    '''
    dtor=np.pi/180.
    min_ele = np.min(rsHH['elepos']/1000.)
    meanele = np.mean(rsHH['elepos']/1000.) 	
    max_ele = np.max(rsHH['elepos']/1000.)
    stadele = np.std(rsHH['elepos']/1000.)
    print(min_ele,max_ele,meanele,stadele)
    rsHH['elepos'][rsHH['elepos']==0]=1000.
    trkmsk=((rsHH['opmode']==0) | (rsHH['opmode']==2) | (rsHH['opmode']==3) | (rsHH['opmode']==5) | (rsHH['opmode']==9))  & ((rsHH['target']/32!=1) & (rsHH['target']/32!=2))
    varele = (rsHH['elepos']/1000.)
    expo212_HH = opa212/np.sin(dtor*varele)
    expo405_HH = opa405/np.sin(dtor*varele)
    aux212_HH = np.exp(-1.0*expo212_HH)
    aux405_HH = np.exp(-1.0*expo405_HH)
    for i  in range(0,4):
        tsky212 = ts212
        #resv1 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term1 = rsHH['adc'][:,i].astype('float')
        term2 = tsky212*(1.-aux212_HH)
        term3 = (term1-term2)/aux212_HH
        #resv3 = np.nanmin(term3)
        #term3 = term3#-resv3+resv1
        #difterm212 = np.where(term1 < term2,term1-0.9*term2,term1-term2)
        #cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 3) | (rsHH['opmode']  == 4) | (rsHH['opmode']  == 5)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        #cond = cond_1 & cond_2
        #tsource212 = np.where(cond,difterm212/aux212_HH,term1)
        #term3.clip(np.nanmin())
        rsHH['adc'][:,i] = term3.astype('int32')
        ##min00 = np.min(rsHH['adcval'][rsHH['opmode']==0,i])
        #if (ceutemp==1):
        #    min10 = np.min(rsHH['adcval'][rsHH['opmode']==10,i])
        #    min99 = np.min(rsHH['adcval'][rsHH['opmode']==99,i])
        #    rsHH['adcval'][:,i] -= np.min([min10,min99])
        #(rsHH['adcval'][:,i]).clip(0.)
    for i  in range(4,6):
        tsky405 = ts405
        #resv4 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term4 = rsHH['adc'][:,i].astype('float')
        term5 = tsky405*(1.-aux405_HH)
        term6 = (term4-term5)/aux405_HH
        #resv6 = np.nanmin(term6)
        #term6 = term6-resv6+resv4
        #difterm405 = np.where(term4 <= term5,term4-0.9*term5,term4-term5)
        #cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 3) | (rsHH['opmode'] == 4) | (rsHH['opmode']  == 5)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        #cond = cond_1 & cond_2
        #tsource405 = np.where(cond,difterm405/aux405_HH,term4)
        #term6.clip()pseudocal
        rsHH['adc'][:,i] = term6.astype('int32')
        ##min00 = np.min(rsHH['adcval'][rsHH['opmode']==0,i])
        #if (ceutemp ==1):
        #    min10 = np.min(rsHH['adcval'][rsHH['opmode']==10,i])
        #    min99 = np.min(rsHH['adcval'][rsHH['opmode']==99,i])
        #    rsHH['adcval'][:,i] -= np.min([min10,min99])
        #(rsHH['adcval'][:,i]).clip(0.)
    return rsHH

def opadescsimple(rsHH,opa212=0.33,ts212=280.,opa405=1.67,ts405=265.0):
    '''
    Descuento de opacidad atmosferica, los tau son ingresados
    del site http://190.3.114.98:20081/sst/sst_plot.php/
    Requiere tambien la tercera taxonomia: iiitaHH (despues
    de la calibracion en temperatura) para modificar los datos
    rsHH, de lo cuales tambien toma las elevaciones. Tiene que
    existir calculo de opacidad en ese archivo rs.
    '''
    dtor=np.pi/180.
    min_ele = np.min(rsHH['elepos']/1000.)
    meanele = np.mean(rsHH['elepos']/1000.) 	
    max_ele = np.max(rsHH['elepos']/1000.)
    stadele = np.std(rsHH['elepos']/1000.)
    print(min_ele,max_ele,meanele,stadele)
    rsHH['elepos'][rsHH['elepos']==0]=1000.
    trkmsk=((rsHH['opmode']==0) | (rsHH['opmode']==2) | (rsHH['opmode']==5) | (rsHH['opmode']==9))  & ((rsHH['target']/32!=1) & (rsHH['target']/32!=2))
    varele = (rsHH['elepos']/1000.)
    expo212_HH = opa212/np.sin(dtor*varele)
    expo405_HH = opa405/np.sin(dtor*varele)
    aux212_HH = np.exp(-1.0*expo212_HH)
    aux405_HH = np.exp(-1.0*expo405_HH)
    for i  in range(0,4):
        tsky212 = ts212
        #resv1 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term1 = rsHH['adcval'][:,i].astype('float')
        term2 = tsky212*(1.-aux212_HH)
        term3 = (term1-term2)/aux212_HH
        #resv3 = np.nanmin(term3)
        #term3 = term3#-resv3+resv1
        #difterm212 = np.where(term1 < term2,term1-0.9*term2,term1-term2)
        #cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 3) | (rsHH['opmode']  == 4) | (rsHH['opmode']  == 5)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        #cond = cond_1 & cond_2
        #tsource212 = np.where(cond,difterm212/aux212_HH,term1)
        #term3.clip(np.nanmin())
        rsHH['adcval'][:,i] = term3.astype('int32')
        ##min00 = np.min(rsHH['adcval'][rsHH['opmode']==0,i])
        #if (ceutemp==1):
        #    min10 = np.min(rsHH['adcval'][rsHH['opmode']==10,i])
        #    min99 = np.min(rsHH['adcval'][rsHH['opmode']==99,i])
        #    rsHH['adcval'][:,i] -= np.min([min10,min99])
        #(rsHH['adcval'][:,i]).clip(0.)
    for i  in range(4,6):
        tsky405 = ts405
        #resv4 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term4 = rsHH['adcval'][:,i].astype('float')
        term5 = tsky405*(1.-aux405_HH)
        term6 = (term4-term5)/aux405_HH
        #resv6 = np.nanmin(term6)
        #term6 = term6-resv6+resv4
        #difterm405 = np.where(term4 <= term5,term4-0.9*term5,term4-term5)
        #cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 3) | (rsHH['opmode'] == 4) | (rsHH['opmode']  == 5)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        #cond = cond_1 & cond_2
        #tsource405 = np.where(cond,difterm405/aux405_HH,term4)
        #term6.clip()pseudocal
        rsHH['adcval'][:,i] = term6.astype('int32')
        ##min00 = np.min(rsHH['adcval'][rsHH['opmode']==0,i])
        #if (ceutemp ==1):
        #    min10 = np.min(rsHH['adcval'][rsHH['opmode']==10,i])
        #    min99 = np.min(rsHH['adcval'][rsHH['opmode']==99,i])
        #    rsHH['adcval'][:,i] -= np.min([min10,min99])
        #(rsHH['adcval'][:,i]).clip(0.)
    return rsHH

def opadesc(rsHH,opa212=0.33,ts212=280.,opa405=1.67,ts405=265.0,deltat212=0.,deltat405=0.,ceutemp=1):
    '''
    Descuento de opacidad atmosferica, los tau son ingresados
    del site http://190.3.114.98:20081/sst/sst_plot.php/
    Requiere tambien la tercera taxonomia: iiitaHH (despues
    de la calibracion en temperatura) para modificar los datos
    rsHH, de lo cuales tambien toma las elevaciones. Tiene que
    existir calculo de opacidad en ese archivo rs.
    '''
    dtor=np.pi/180.
    min_ele = np.min(rsHH['elepos']/1000.)
    max_ele = np.max(rsHH['elepos']/1000.)
    stadele = np.std(rsHH['elepos']/1000.)
    print (min_ele,max_ele,stadele)
    ###############################################
    scantauori=rsHH['adcval'][rsHH['opmode']==10,:]
    scanintori=rsHH['adcval'][rsHH['opmode']==9,:]
    ############################################### 
    rsHH['elepos'][rsHH['elepos']==0]=1000.
    trkmsk=((rsHH['opmode']==0) | (rsHH['opmode']==2) | (rsHH['opmode']==3) | (rsHH['opmode']==4) | (rsHH['opmode']==5) | (rsHH['opmode']==9) | (rsHH['opmode']==50) | (rsHH['opmode']==99))  & ((rsHH['target']/32!=1) & (rsHH['target']/32!=2))
    varele = (rsHH['elepos'][trkmsk]/1000.)
    expo212_HH = opa212/np.sin(dtor*varele)
    expo405_HH = opa405/np.sin(dtor*varele)
    aux212_HH = np.exp(-1.0*expo212_HH)
    aux405_HH = np.exp(-1.0*expo405_HH)
    for i  in range(0,4):
        tsky212 = ts212
        resv1 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term1 = rsHH['adcval'][trkmsk,i].astype('float')
        term2 = tsky212*(1.-aux212_HH)
        term3 = (term1-term2)/aux212_HH
        resv3 = np.nanmin(term3)
        term3 = term3-resv3
        rsHH['adcval'][trkmsk,i] = term3.astype('int32')
        rsHH['adcval'][rsHH['opmode']==10,i] = scantauori[:,i]
        rsHH['adcval'][rsHH['opmode']== 9,i] = scanintori[:,i]  
    for i  in range(4,6):
        tsky405 = ts405
        resv4 = np.nanmin(rsHH['adcval'][trkmsk,i].astype('float'))
        term4 = rsHH['adcval'][trkmsk,i].astype('float')
        term5 = tsky405*(1.-aux405_HH)
        term6 = (term4-term5)/aux405_HH
        resv6 = np.nanmin(term6)
        term6 = term6-resv6
        rsHH['adcval'][trkmsk,i] = term6.astype('int32')
        rsHH['adcval'][rsHH['opmode']==10,i] = scantauori[:,i]
        rsHH['adcval'][rsHH['opmode']== 9,i] = scanintori[:,i] 
    return rsHH

def efiefe(ivtaHH,rsHH,corr=0):
    '''
    Calcula la eficiencia efectiva por haz($\eta_{i}$) segun
    la ecuacion (1) de Krucker et al 2013, asumiendo tempera-
    turas de Sol calmo de 5900 K y 5100 K para 212 y 405 res-
    pectivamente segundo Silva et al 2005.
    '''
    trkmsk=((rsHH['opmode']==0) | (rsHH['opmode']==2) | (rsHH['opmode']==5) | (rsHH['opmode']==9)) & ((rsHH['target']/32!=1) & (rsHH['target']/32!=2))
    print( len(trkmsk))
    beff = np.zeros(shape=(6,1),dtype=float)
    for i  in range(0,4):
        infe212 = np.nanmin(ivtaHH['scnval'][:,i])
        supe212scn = np.nanmax(ivtaHH['scnval'][:,i])
        supe212trk = np.nanmax(ivtaHH['trkval'][:,i])
        if corr==1:
            supe212 = supe212scn
        elif corr==0:
            supe212 = supe212trk 
        else:
            print("Are you in drugs?")
        efic212 = (supe212-infe212)/(5900.)
        beff[i] = efic212
        print(efic212)
        rsHH['adcval'][trkmsk,i] -= infe212
        rsHH['adcval'][trkmsk,i] = 5900. + rsHH['adcval'][trkmsk,i]/efic212 
    for i in range(4,6):
        infe405 = np.nanmin(ivtaHH['scnval'][:,i])
        supe405scn = np.nanmax(ivtaHH['scnval'][:,i])
        supe405trk = np.nanmax(ivtaHH['trkval'][:,i])
        if corr==1:
            supe405 = supe405scn
        elif corr==0:
            supe405 = supe405trk
        else:
            print("No option")
        efic405 = (supe405-infe405)/(5100.)
        beff[i] = efic405
        print(efic405)
        rsHH['adcval'][trkmsk,i] -= infe405
        rsHH['adcval'][trkmsk,i] = 5100. + rsHH['adcval'][trkmsk,i]/efic405 
    print('Quiet-Sun brightness temperature 5900K@212GHz,5100K@405GHz, Silva et al (2005)')
    return rsHH,beff

def beameffmap(ivtaHH,rsHH,corr=0):
    '''
    Calcula la eficiencia efectiva por haz($\eta_{i}$) segun
    la ecuacion (1) de Krucker et al 2013, asumiendo tempera-
    turas de Sol calmo de 6344 K y 5858 K para 212 y 405 res-
    pectivamente.
    '''
    beff = np.zeros(shape=(6,1),dtype=float)
    for i  in range(0,4):
        #infe212scn = np.min(ivtaHH['scnval'][:,i])
        #infe212tau = np.min(ivtaHH['tauval'][:,i])
        infe212 = np.min(ivtaHH['mapval'][:,i])#np.where(infe212scn <= infe212tau,infe212scn,infe212tau)
        supe212scn = np.max(ivtaHH['mapval'][:,i])
        supe212trk = np.max(ivtaHH['trkval'][:,i])
        if corr==1:
            supe212 = np.where(supe212scn <= supe212trk,supe212scn,supe212trk)
        else:
            supe212 = supe212trk #np.where(supe212scn >= supe212trk,supe212scn,supe212trk)
        efic212 = (supe212)/(6344.)
        beff[i] = efic212
        print (efic212)
        cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 5) | (rsHH['adcval'][:,i] >= supe212/2.)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        cond = cond_1 & cond_2
        rsHH['adcval'][:,i] = np.where(cond, rsHH['adcval'][:,i]/efic212, rsHH['adcval'][:,i])
    for i in range(4,6):
        #infe405scn = np.min(ivtaHH['scnval'][:,i])
        #infe405tau = np.min(ivtaHH['tauval'][:,i])
        infe405 =  np.min(ivtaHH['mapval'][:,i])#np.where(infe405scn <= infe405tau,infe405scn,infe405tau)
        supe405scn = np.max(ivtaHH['mapval'][:,i])
        supe405trk = np.max(ivtaHH['trkval'][:,i])
        if corr==1:
            supe405 = np.where(supe405scn <= supe405trk,supe405scn,supe405trk)
        else:
            supe405 = supe405trk
        efic405 = (supe405)/(5858.)
        beff[i] = efic405
        print(efic405)
        cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 5) | (rsHH['adcval'][:,i] >= supe405/2.)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        cond = cond_1 & cond_2
        rsHH['adcval'][:,i] = np.where(cond, rsHH['adcval'][:,i]/efic405, rsHH['adcval'][:,i])
    for i in range(0,4):
        rsHH['adcval'][:,i]=np.mean(beff[0:4])*rsHH['adcval'][:,i]
    for i in range(4,6):
        rsHH['adcval'][:,i]=np.mean(beff[4:6])*rsHH['adcval'][:,i]
    return rsHH,beff

def beameffrdm(ivtaHH,rsHH,corr=0):
    '''
    Calcula la eficiencia efectiva por haz($\eta_{i}$) segun
    la ecuacion (1) de Krucker et al 2013, asumiendo tempera-
    turas de Sol calmo de 6344 K y 5858 K para 212 y 405 res-
    pectivamente.
    '''
    beff = np.zeros(shape=(6,1),dtype=float)
    for i  in range(0,4):
        #infe212scn = np.min(ivtaHH['scnval'][:,i])
        #infe212tau = np.min(ivtaHH['tauval'][:,i])
        infe212 = np.min(ivtaHH['rdmval'][:,i])#np.where(infe212scn <= infe212tau,infe212scn,infe212tau)
        supe212scn = np.max(ivtaHH['rdmval'][:,i])
        supe212trk = np.max(ivtaHH['rdmval'][:,i])
        if corr==1:
            supe212 = np.where(supe212scn <= supe212trk,supe212scn,supe212trk)
        else:
            supe212 = supe212trk #np.where(supe212scn >= supe212trk,supe212scn,supe212trk)
        efic212 = (supe212)/(6344.)
        beff[i] = efic212
        print(efic212)
        cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 5) | (rsHH['adcval'][:,i] >= supe212/2.)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        cond = cond_1 & cond_2
        rsHH['adcval'][:,i] = np.where(cond, rsHH['adcval'][:,i]/efic212, rsHH['adcval'][:,i])
    for i in range(4,6):
        #infe405scn = np.min(ivtaHH['scnval'][:,i])
        #infe405tau = np.min(ivtaHH['tauval'][:,i])
        infe405 =  np.min(ivtaHH['rdmval'][:,i])#np.where(infe405scn <= infe405tau,infe405scn,infe405tau)
        supe405scn = np.max(ivtaHH['rdmval'][:,i])
        supe405trk = np.max(ivtaHH['rdmval'][:,i])
        if corr==1:
            supe405 = np.where(supe405scn <= supe405trk,supe405scn,supe405trk)
        else:
            supe405 = supe405trk
        efic405 = (supe405)/(5858.)
        beff[i] = efic405
        print (efic405)
        cond_1 = (rsHH['opmode'] == 0) | (rsHH['opmode'] == 2) | (rsHH['opmode'] == 5) | (rsHH['adcval'][:,i] >= supe405/2.)
        #cond_2 = (rsHH['target']/32 <> 1) & (rsHH['target']/32 <> 2)
        cond = cond_1 & cond_2
        rsHH['adcval'][:,i] = np.where(cond, rsHH['adcval'][:,i]/efic405, rsHH['adcval'][:,i])
    for i in range(0,4):
        rsHH['adcval'][:,i]=np.mean(beff[0:4])*rsHH['adcval'][:,i]
    for i in range(4,6):
        rsHH['adcval'][:,i]=np.mean(beff[4:6])*rsHH['adcval'][:,i]
    return rsHH,beff

def convol(Xplusone,beamX,veclen,pos,mr,fi,tinf=273.,tsup=6180.,delta=0,smu=1):
    '''

    '''
    X = Xplusone-1
    a = campo(0.,1.,mr)
    xcalX = (pos.beams_ew[fi,X]+0.)/3.6+len(a)/2
    ycalX = (pos.beams_ns[fi,X]+0.)/3.6+len(a)/2
    convsumX = np.zeros(len(xcalX))
    for j in range(0,veclen):
        unverX = a[ycalX[j].astype(int)-300+delta:ycalX[j].astype(int)+301-delta:1,xcalX[j].astype(int)-300+delta:xcalX[j].astype(int)+301-delta:1]
        convsumX[j] = np.sum(unverX*beamX[delta:601-delta,delta:601-delta])
    df  = pd.DataFrame(data=convsumX)
    convsumX = df.rolling(min_periods=1,center=True,window=smu).mean()
    convsumX = (convsumX-1.*np.min(convsumX))/np.max(convsumX)
    convsumX = tsup*convsumX+tinf
    #convsumX -= np.nanmin(convsumX)
    return convsumX

def convoltest(Xplusone,beamX,veclen,pos,mr,fi,t_ini=273.,t_fin=5900.,smu=1,delta=0):
    '''
    
    '''
    X = Xplusone-1
    a = campo(0.,1.,mr)
    xcalX = (pos.beams_ew[fi,X]+0.)/3.6+len(a)/2
    ycalX = (pos.beams_ns[fi,X]+0.)/3.6+len(a)/2
    convsumX = np.zeros(len(xcalX))
    beamX = beamX/beamX.sum()
    for j in range(0,veclen):
        #a[arcy-arra:arcy+arra, arcx-arra:arcx+arra] = ar[:,:,i]
        unverX = a[ycalX[j].astype(int)-300+delta:ycalX[j].astype(int)+301-delta:1,\
                   xcalX[j].astype(int)-300+delta:xcalX[j].astype(int)+301-delta:1]
        convsumX[j] = np.sum(unverX*beamX[delta:601-delta,delta:601-delta])
    convsumXdf = pd.DataFrame(convsumX)###,smu,center=False,min_periods=0)
    aux = convsumXdf.rolling(smu,center=True,min_periods=0).sum()
    nconvsumX = (aux[0]-np.min(aux[0]))/np.max(aux[0])
    #convsumX[0] = t_ini
    nconvsumX = (t_fin-t_ini)*nconvsumX+t_ini
    #convsumX -= np.nanmin(convsumX) 
    return nconvsumX


def oldconvol(Xplusone,beamX,veclen,pos,mr,fi,tinf=273,tsup=6340,smu=1,delta=0):
    '''

    '''
    X = Xplusone-1
    a = campo(0.,1.,mr)
    xcalX = (2*pos.source_ew[fi]-pos.beams_ew[fi,X]+0.)/3.6+len(a)/2
    ycalX = (2*pos.source_ns[fi]-pos.beams_ns[fi,X]+0.)/3.6+len(a)/2
    convsumX = np.zeros(len(xcalX))
    for j in range(0,veclen):
        unverX = a[ycalX[j].astype(int)-300+delta:ycalX[j].astype(int)+301-delta:1,xcalX[j].astype(int)-300+delta:xcalX[j].astype(int)+301-delta:1]
        convsumX[j] = np.sum(unverX*beamX[delta:601-delta,delta:601-delta])
    df = pd.DataFrame(data=convsumX)
    convsumX = df.rolling(min_periods=1,center=True,window=smu).mean()
    convsumX = (convsumX-1.0*np.min(convsumX))/np.max(convsumX)
    convsumX = tsup*convsumX+tinf
    return convsumX


#2*pos.source_ew[fi]-pos.beams_ew[fi,j],2*pos.source_ns[fi]-pos.beams_ns[fi,j]


def convollimbo(Xplusone,beamX,veclen,pos,mr,fi,tceu=273,tsol=6340,smu=1,delta=0):
    '''degree

    '''
    X = Xplusone-1
    a = campo(0.,1.,mr)
    xcalX = (pos.beams_ew[fi,X]+10.)/3.6+len(a)/2
    ycalX = (pos.beams_ns[fi,X]+0.)/3.6+len(a)/2
    convsumX = np.zeros(len(xcalX))
    for j in range(0,veclen):
        unverX = a[ycalX[j].astype(int)-300+delta:ycalX[j].astype(int)+301-delta:1,xcalX[j].astype(int)-300+delta:xcalX[j].astype(int)+301-delta:1]
        convsumX[j] = np.sum(unverX*beamX[delta:601-delta,delta:601-delta])
    #convsumX = pd.rolling_mean(convsumX,smu,center=True,min_periods=0)
    convsumX = tsol*convsumX/np.max(convsumX)+tceu
    return convsumX

def convolmap(veclen,a,xcalI,ycalI,beamI,delta=0):
    '''

    '''
    convsum = np.zeros(veclen)
    for j in range(0,veclen):
        unver = a[ycalI[j].astype(int)-300+delta:ycalI[j].astype(int)+301-delta,\
                  xcalI[j].astype(int)-300+delta:xcalI[j].astype(int)+301-delta]
        pizza = beamI[delta:601-delta,delta:601-delta]
        convsum[j] = np.sum(unver*pizza)
    return convsum

def convolmaptotal(veclen,a,xcalI,ycalI,beamI,delta=0):
    '''

    '''
    convsum = np.zeros(veclen)
    for j in range(0,veclen):
        unver = a[ycalI[j]-300+delta:ycalI[j]+301-delta,xcalI[j]-300+delta:xcalI[j]+301-delta]
        pizza = beamI[delta:601-delta,delta:601-delta]
        convsum[j] = np.sum(unver*pizza)
    return convsum

def diffcalc(Xplusone,rsHH,iniz,fina,convsumX,alf=True,cut=1.0,adcorr=1.0):
    '''
    Ajuste lineal de las diferencias entre observacion y
    convolucion, considerando apenas la parte de "tracking"
    La correcion tiene efecto si la variable alf==True
    '''
    X = Xplusone-1
    ffina = np.int(cut*fina)
    obsadcX = rsHH['adcval'][:,X][iniz:fina]# - np.min(rsHH['adcval'][:,0][iniz:fina]).astype(float)
    normcs0X = np.mean(obsadcX[rsHH['opmode']==0][iniz:ffina])*convsumX/np.mean(convsumX[rsHH['opmode']==0][iniz:ffina])
    d0X = (obsadcX-normcs0X)
    residual_std_error = 0.
    ###################################################################################################
    # Ajuste lineal de las diferencias usando scipy.stats
    if alf==True:
        slope,intercept,r_value,p_value,std_err=stats.linregress((rsHH['opmode'][iniz:ffina]==0)*(rsHH['time'][iniz:ffina]),d0X[iniz:ffina])
        #Calculate some additional outputs
        predict_y = intercept + slope * (rsHH['opmode'][iniz:ffina]==0)*(rsHH['time'][iniz:ffina])
        pred_error = d0X[iniz:ffina] - predict_y
        degrees_of_freedom = len((rsHH['opmode'][iniz:ffina]==0)*(rsHH['time'][iniz:ffina])) - 2
        residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
        normcsX = normcs0X + intercept*(rsHH['opmode']==0) + slope*(rsHH['opmode']==0)*(rsHH['time'][iniz:fina])
        dX = (obsadcX-normcsX)
    else:
        normcsX = normcs0X
        dX = d0X
    sumX = np.sqrt(np.sum(dX**2)/np.sum(obsadcX**2))
    dX = adcorr*dX # force an ad libitum correction
    ###################################################################################################
    return dX,sumX,residual_std_error

def detrendme(rsHH,mode=0,ini=0,fin=-1,channels=[0,1,2,3,4,5],window=25,dt=True):
    '''
    Correccion de los canales de rs contenidos en la var. channels
    aplicando un "detrend" linear. Para descartar la correccion y
    solo plotar hacer dt=False.
    '''
    #import matplotlib
    #bpw=[a/window for a in BP]
    indice = np.linspace(0,len(rsHH['time'])-1,len(rsHH['time']))
    #shortindex=np.where(rscc['adcval']==rscc['adcval'][np.where(rscc['opmode']==10)][0])
    x0,x1=indice[ini],indice[fin]
    y0,y1=rsHH['adcval'][ini],rsHH['adcval'][fin]
    for i in channels:
        z = np.polyfit([x0,x1],[y0[i],y1[i]], 1)
        p = np.poly1d(z)
        aux = rsHH['adcval'][:,i] - p(indice).astype('int32')
        aux -= np.min(aux)
        if dt==False:
            plt.plot(aux)
        else:
           rsHH['adcval'][:,i]=aux
           #plt.plot(aux)
        #plt.ylim(12000,25000)
        #plt.grid()
    return rsHH

def cutme(tihh,rsHH,poshh,w,cpi,cpf):
    '''


    '''
    INI,FIN=np.floor(cpi/w),np.floor(cpf/w)
    for items in rsHH:
        rsHH[items]=rsHH[items][INI:FIN]
    tihh = tihh[INI:FIN]
    for items in poshh:
        poshh[items]=poshh[items][INI:FIN]
    return tihh,rsHH,poshh


##########################################################################################################
##########################################################################################################
################################################MULTIBEAM#################################################

def multibeam(rs,equis=[0,0,0,0,0,0],ygrie=[0,0,0,0,0,0],ef212=0.20,ef405=0.10,tqs=5900.):
    '''



    '''
    #POSIcaO DOS FEIXES
    xideal=[-0.5,1.5,-2.5,-0.5,0.0,-0.5]
    yideal=[7.5,1.5,1.5,-2.0,0.0,7.5]

    x0=13.9799999*60. + equis[0]
    x1=16.0499992*60. + equis[1]
    x2=12.1199999*60. + equis[2]
    x3=14.2500000*60. + equis[3]
    x4=14.5299997*60. + equis[4]
    x5=13.8599997*60. + equis[5]

    y0=-4.44999980*-60. + ygrie[0]
    y1=-10.6199999*-60. + ygrie[1]
    y2=-10.4600000*-60. + ygrie[2]
    y3=-14.3500004*-60. + ygrie[3]
    y4=-12.1099997*-60. + ygrie[4]
    y5=-4.34000020*-60. + ygrie[5]

    x = [x0,x1,x2,x3,x4,x5]
    y = [y0,y1,y2,y3,y4,y5]


    #DESVIOS PADRaO
    sigma212 = 1.6*60.*4.1/np.sqrt(np.log(256)) # 212 GHz
    sigma405_a = 60.*2.9/np.sqrt(np.log(256)) # 405 GHz
    sigma405_b = 60.*5.2/np.sqrt(np.log(256)) # 

    sigma212_2 = sigma212*sigma212  # square(212)
    sigma405_2 = sigma405_a*sigma405_b  # square(405)

    sigma2x212_2 = 2.*sigma212_2 # 2*square(212)
    sigma2x405_2 = 2.*sigma405_2 # 2*square(405)

    numer = np.size(rs['adcval'][:,0]) # la longitud del trecho
    T_Ant = rs['adcval'].astype(float) # seis canales, 4@212 y 2@405
    T_Ant[:,1] = tqs + T_Ant[:,1] - np.nanmin(T_Ant[:,1])
    T_Ant[:,2] = tqs + T_Ant[:,2] - np.nanmin(T_Ant[:,2])
    T_Ant[:,3] = tqs + T_Ant[:,3] - np.nanmin(T_Ant[:,3])

    Tca = np.ones(shape=(numer,4),dtype=float) #cuatro canales 3@212 y 1@405

    #Logaritmo del ratio de temperaturas
    lr_T23 = np.log(T_Ant[:,1]/T_Ant[:,2]) #log de (T_Chan2/T_Chan3)
    lr_T34 = np.log(T_Ant[:,2]/T_Ant[:,3]) #log de (T_Chan3/T_Chan4)
    lr_T42 = np.log(T_Ant[:,3]/T_Ant[:,1]) #log de (T_Chan4/T_Chan2)

    #Parametro de contraste K
    K=lr_T42-lr_T23 #K=ln(TH*TI/TL^2)

    #CALCULO DE multiplos feixes
    #SOMA E RESTA DE COORDENADAS
    sx12=x1+x2
    rx12=x1-x2
    sx13=x1+x3
    rx13=x1-x3
    sy12=y1+y2
    ry12=y1-y2
    sy13=y1+y3
    ry13=y1-y3
    sx23=x2+x3
    rx23=x2-x3
    sy23=y2+y3
    ry23=y2-y3

    #numerador de phi
    A=-2.*rx12*sigma212_2*lr_T42+rx13*sx13*rx12+ry13*sy13*rx12-2.*rx13*sigma212_2*lr_T23-rx12*sx12*rx13-ry12*sy12*rx13

    #denominador de phi
    B=2.*ry13*rx12-2.*ry12*rx13

    #phi
    sy=A/B

    #numerador de theta
    C=2.*sigma212_2*lr_T23+rx12*sx12+ry12*sy12-2.*sy*ry12

    #denominador de theta
    D=2.*rx12

    #theta
    sx=C/D

    #Reconstruccion da temperatura/
    for i in range(1,4):
        diff212x2=(x[i]-sx)*(x[i]-sx)
        diff212y2=(y[i]-sy)*(y[i]-sy)
        Tca[:,i-1] = tqs + (T_Ant[:,i] - np.nanmin(T_Ant[:,i]))*np.exp((diff212x2+diff212y2)/sigma2x212_2)

    #Tentativa de reconstruccion de Tb @ 405 Ghz
    diff405x2=(x[4]-sx)*(x[4]-sx)
    diff405y2=(y[4]-sy)*(y[4]-sy)
    Tca[:,3] = tqs*5.1/5.9 + (T_Ant[:,4] - np.nanmin(T_Ant[:,4]))*np.exp((diff405x2+diff405y2)/sigma2x405_2)


    #Calculo de los fluxos
    Flux = np.ones(shape=(numer,2),dtype=float)
    kb = 1.38e-23
    Ageom = 1.5*1.5*np.pi/4.
    Aeff212 = ef212*Ageom
    Aeff405 = ef405*Ageom
    Flux[:,0] = 2.*kb*Tca[:,1]/Aeff212*1.e22 #212
    Flux[:,0] = Flux[:,0]-np.min(Flux[:,0])
    Flux[:,1] = 2.*kb*Tca[:,3]/Aeff405*1.e22 #405
    Flux[:,1] = Flux[:,1]-np.min(Flux[:,1])
    return Flux,Tca,sx-x4,sy-y4,K

##########################################################################################################
##########################################################################################################
##########################################################################################################


##########################################################################################################
##########################################################################################################
################################################MULTIBEAM#################################################

def multibeamtest(rs,pos,ef212=0.35,ef405=0.17,errotemp=20.,gaussf=1.1026):
    '''



    '''
    #POSIcaO DOS FEIXES


    x0=pos.beams_ew[:,0]
    x1=pos.beams_ew[:,1]
    x2=pos.beams_ew[:,2]
    x3=pos.beams_ew[:,3]
    x4=pos.beams_ew[:,4]
    x5=pos.beams_ew[:,5]

    y0=pos.beams_ns[:,0]
    y1=pos.beams_ns[:,1]
    y2=pos.beams_ns[:,2]
    y3=pos.beams_ns[:,3]
    y4=pos.beams_ns[:,4]
    y5=pos.beams_ns[:,5]

    x = [x0,x1,x2,x3,x4,x5]
    y = [y0,y1,y2,y3,y4,y5]


    #DESVIOS PADRaO
    sigma212 = 60.*3.9750*gaussf/np.sqrt(np.log(256.)) # 212 GHz
    #sigma212_b = 60.*4.0625*gaussf/np.sqrt(np.log(256.))
    sigma405 = 60.*2.9/np.sqrt(np.log(256.)) # 405 GHz
    #sigma405_b = 60.*5.2/np.sqrt(np.log(256.)) # 

    sigma212_2 = sigma212*sigma212  # square(212)
    #sigma212_2b = sigma212_b*sigma212_b
    sigma405_2 = sigma405*sigma405  # square(405)
    #sigma405_2b = sigma405_b*sigma405_b  

    sigma2x212_2 = 2.*sigma212_2 # 2*square(212)
    #sigma2x212_2b = 2.*sigma212_2b
    sigma2x405_2 = 2.*sigma405_2 # 2*square(405)
    #sigma2x405_2b = 2.*sigma405_2b

    numer = np.size(rs['adcval'][:,0]) # la longitud del trecho
    T_Ant = rs['adcval'].astype(float) # seis canales, 4@212 y 2@405
    Tbkg212 = T_Ant[:,0]
    Tbkg405 = T_Ant[:,5]
    T_Ant[:,1] = np.abs(T_Ant[:,1] - T_Ant[:,0]+errotemp)
    T_Ant[:,2] = np.abs(T_Ant[:,2] - T_Ant[:,0]+errotemp)
    T_Ant[:,3] = np.abs(T_Ant[:,3] - T_Ant[:,0]+errotemp)
    T_Ant[:,4] = np.abs(T_Ant[:,4] - T_Ant[:,5]+errotemp)

    Tca = np.ones(shape=(numer,4),dtype=float) #cuatro canales 3@212 y 1@405

    #Logaritmo del ratio de temperaturas
    lr_T23 = np.log(T_Ant[:,1]/T_Ant[:,2]) #log de (T_Chan2/T_Chan3)
    lr_T34 = np.log(T_Ant[:,2]/T_Ant[:,3]) #log de (T_Chan3/T_Chan4)
    lr_T42 = np.log(T_Ant[:,3]/T_Ant[:,1]) #log de (T_Chan4/T_Chan2)

    #Parametro de contraste K
    K=lr_T42-lr_T23 #K=ln(TH*TI/TL^2)

    #CALCULO DE multiplos feixes
    #SOMA E RESTA DE COORDENADAS
    sx12=x1+x2
    rx12=x1-x2
    sx13=x1+x3
    rx13=x1-x3
    sy12=y1+y2
    ry12=y1-y2
    sy13=y1+y3
    ry13=y1-y3
    sx23=x2+x3
    rx23=x2-x3
    sy23=y2+y3
    ry23=y2-y3
#'''
    #numerador de phi
    A = -2. * rx12 * sigma212_2 * lr_T42 + 1. * rx13 * sx13 * rx12 + 1. * ry13 * sy13 * rx12 -2. * rx13 * sigma212_2 * lr_T23 - rx12 * sx12 * rx13 - 1. * ry12 * sy12 * rx13

    #denominador de phi
    B = 2. * ry13 * rx12 - 2. * ry12 * rx13
 
    #phi
    sy=(A/B)
   
    #numerador de theta
    C= 2. * sigma212_2 * lr_T23 + 1. * rx12 * sx12+1. * ry12 * sy12-1. * sy * ry12

    #denominador de theta
    D= 2. * rx12
   
    #theta
    sx=(C/D)
  
    #Reconstruccion da temperatura/
    for i in range(1,4):
        diff212x2 = (sx-x[i])**2
        diff212y2 = (sy-y[i])**2
        Tca[:,i-1] = T_Ant[:,i]*np.exp((diff212x2+diff212y2)/sigma2x212_2)

    #Tentativa de reconstruccion de Tb @ 405 Ghz
    diff405x2=(sx-x[4])**2
    diff405y2=(sy-y[4])**2
    Tca[:,3] = T_Ant[:,4]*np.exp((diff405x2+diff405y2)/sigma2x405_2)


    #Calculo de los fluxos
    Flux = np.ones(shape=(numer,2),dtype=float)
    kb = 1.38e-23
    Ageom = 1.5*1.5*np.pi/4.
    Aeff212 = ef212*Ageom
    Aeff405 = ef405*Ageom
    Flux[:,0] = 2.*kb*Tca[:,1]/Aeff212*1.e22 #212
    #Flux[:,0] = Flux[:,0]-np.min(Flux[:,0])
    Flux[:,1] = 2.*kb*Tca[:,3]/Aeff405*1.e22 #405
    #Flux[:,1] = Flux[:,1]-np.min(Flux[:,1])
    return Tca,sx,sy,K,Tbkg212,Tbkg405
#'''
##########################################################################################################
##########################################################################################################
##########################################################################################################

def  aarcalc(rs,rel,semilado,visrad):
    '''
    Alt-Azimuth maps radius calculation
    '''

    for nume in range(0,6):
        x2 = rs['opmode'] == 2 # 3 es el tag que describe el mapa radial
        newx2 = ((rs['time'] >= np.min(rs['time'][x2])) & (rs['time'] <= np.max(rs['time'][x2])))
        #dif2x2newx2 = np.where(x2 <> newx2)
        x2 = newx2
        xoff,yoff = rs['x_off'][x2],rs['y_off'][x2] # de los datos originales, que luegos seran
        veclen = len(xoff) # longitud/tamanho del vector posicion
        iniz = 0
        fina = veclen
        xbeampos,ybeampos = rel.bpos.off[0][nume]*60./3.6,rel.bpos.el[0][nume]*60./3.6
        xcenter,ycenter = rel.bpos.off[0][4]*60./3.6,rel.bpos.el[0][4]*60./3.6
        xbeampos -= xcenter # beam 5 como centro
        ybeampos -= ycenter # beam 5 como centro
        xcenter -= xcenter
        ycenter -= ycenter
        xcenlist,ycenlist = xoff+semilado,yoff+semilado # recorrido del centro de apont. desplazado al cuadrante pos.
        xlist,ylist = xcenlist + xbeampos, ycenlist + ybeampos # recorrido del haz 1 desplazado al cuadrante pos.
        ################################################################################
        '''
        plt.figure(figsize=(6,6))
        plt.plot(xbeampos,ybeampos,'*',markersize=23,label='haz')
        plt.plot(xcenter,ycenter,'*',markersize=23,label='centro')
        plt.plot(xcenlist,ycenlist,label='recorr. central',lw=1)
        plt.plot(xlist,ylist,label ='recorr. haz',lw=1)
        plt.text(500,750,'dist.(arc-min) = '+np.str( np.sqrt((xbeampos*3.6)**2+(ybeampos*3.6)**2)/60.))
        plt.text(500,450,'X (arc-min) : '+np.str(xbeampos*3.6/60.))
        plt.text(500,150,'Y (arc-min) : '+np.str(ybeampos*3.6/60.))
        plt.legend(loc=2,frameon=False)
        #plt.xlim(1900,2100)
        #plt.ylim(1900,2100)
        plt.grid()
        '''
        ###########################################################
        delta= 0.05
        inflim,suplim=0.5-delta,0.5+delta
        #############################################################
        obsadc = (rs['adcval'][x2,nume][iniz:fina] - np.min(rs['adcval'][x2,nume][iniz:fina])).astype(float)
        eje = np.linspace(0,len(obsadc)-1,len(obsadc),dtype='int')
        doublemed = np.percentile(obsadc,99.99) ## valores con percentil superior
        fifty = np.where((obsadc <= suplim*doublemed) & (obsadc >=inflim*doublemed))
        #np.size(fifty)
        topbotdif = np.max(obsadc)-np.min(obsadc)
        #primeira determinacao
        [firstmaxtab, firstmintab]=peakdet(obsadc,0.75*topbotdif,eje)
        firstone1 = firstmintab[:,0].astype(int)
        firstone2 = np.array(np.where((obsadc <= 0.0025*doublemed) & (obsadc >=0.0*doublemed)))
        firstonedy = np.concatenate((firstone1,firstone2[0]),axis=0)
        firsthundry = firstmaxtab[:,0].astype(int)
        c = np.polyfit(eje[firstonedy],obsadc[firstonedy],1) #ajuste polinomial de valores base (minimos)
        f = np.poly1d(c)
        oldbase = f(eje)
        nobsadc = obsadc - oldbase + obsadc[firstonedy].mean()
        oldadc= obsadc
        obsadc = nobsadc# - np.min(nobsadc1)
        #segunda determinacao
        [secondmaxtab, secondmintab]=peakdet(obsadc,0.75*topbotdif,eje)
        secondone1 = secondmintab[:,0].astype(int)
        secondone2 = np.where((obsadc <= 0.0025*doublemed) & (obsadc >=0.0*doublemed))
        firstonedy = np.concatenate((secondone1,secondone2[0]),axis=0)
        secondhundry = firstmaxtab[:,0].astype(int) #conjunto valores maximos y ...
        newtwomed = np.percentile(obsadc,99.99)
        nfifty = np.where((obsadc <= suplim*newtwomed) & (obsadc >=inflim*newtwomed))
        #np.size(nfifty)
        '''
        plt.figure(figsize=(11.69*0.5,8.27*.5))
        plt.plot(eje,oldadc,'--',alpha=0.935,label='old')
        plt.plot(eje[firstonedy],oldadc[firstonedy],'*r',ms=18)
        plt.plot(eje[firsthundry],obsadc[firsthundry],'*b',ms=18)
        plt.plot(eje,nobsadc,marker=',',color='orange',label='new')
        plt.plot(eje,oldbase,'g')
        plt.plot(eje[nfifty],nobsadc[nfifty],'xr',ms=13)
        #plt.ylim(-100,500)
        plt.legend(frameon=True)
        #plt.xlim(19577,19599)
        '''
        trasxbeampos, trasybeampos, R2real, residu2 = leastsq_circle(xlist[nfifty],ylist[nfifty])
        myresidualindex = residu2/(np.size(nfifty)-2)
        R2realsec = R2real*3.6
        RelR2 = R2realsec/visrad
        '''
        plot_data_circle(xlist[nfifty],ylist[nfifty],trasxbeampos, trasybeampos, R2real)
        '''
        DeltaX,DeltaY = 3.6*(trasxbeampos-semilado),3.6*(trasybeampos-semilado)
        print('|',nume+1,'|',"%.2f"%DeltaX,'|',"%.2f"%DeltaY,'|',"%.2f"%myresidualindex,'|',"%.2f"%R2realsec,'|',"%.3f"%RelR2,'|')
    return

def J(freq_ghz,T_K,print_error=False,h=h,k_B=k_B):
    '''
    Effective source radiation temperature calculation using Mangum 1993 and Kutner and Ulich 1981
    '''
    freq = freq_ghz*1e9
    nume = (h.value*freq/k_B.value)
    deno = np.exp((h.value*freq)/(k_B.value*T_K)) - 1.
    esrt = nume/deno
    percent_error = 100*np.abs(T_K-esrt)/T_K
    if print_error:
        print('%3.1f'%percent_error)
    return esrt

def etamb_w_moon(freqghz,Ta_star,T_K,theta_mb):
    '''
    Main beam efficiency calculation using Mangum 1993
    '''
    theta_eq,theta_pol = 30.,30.
    polemic = (1 -np.exp(-1*np.log(2.0)*(theta_eq*theta_pol/theta_mb**2)))**-1
    nume_eta_mb = 0.5*(Ta_star)
    deno_eta_mb = J(freqghz,T_K) -J(freqghz,2.79)
    etamb = (nume_eta_mb/deno_eta_mb)*polemic
    return etamb

