#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm ##mapa de colores
from matplotlib import dates

import numpy as np
import pylab
import pandas as pd
import random
from datetime import datetime, timedelta

from scipy.io.idl import readsav
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
from sunpy import lightcurve as lc
from sunpy.time import TimeRange
#from sunpy.net import hek

def t_plot(fds,bi,marker=',',linestyle='solid'):
    """
    plota todos los valores 'adc' de la estrutura 'bi'
    contra el vector-tiempo contenido en 'fds'
    """
    hfmt = dates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(11.69*0.5-8.27*0.5))
    ax = fig.add_subplot(1,1,1)
    b1,b2,b3,b4,b5,b6 = ax.plot(fds,bi.adc,marker=marker,linestyle=linestyle)
    ax.xaxis.set_major_formatter(hfmt)
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_minor_locator(dates.MinuteLocator([15,30,45]))
    ax.xaxis.grid(True,'major',lw=2)
    ax.xaxis.grid(True,'minor',lw=0.5)
    plt.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5","beam 6"], loc='best')
    plt.xticks(rotation=0)
    plt.show()
    return

def bi_plot(fds,bi,marker=',',linestyle='solid'):
    """
    plota todos los valores 'adc' de la estrutura 'bi'
    contra el vector-tiempo contenido en 'fds'
    """
    hfmt = dates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(8.27,8.27/1.618),num=None, dpi=100, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi['adc'][:,0:6],marker=marker,linestyle=linestyle)
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_minor_locator(dates.MinuteLocator([20,40]))
    ax.xaxis.grid(True,'major',lw=2)
    ax.xaxis.grid(True,'minor',lw=0.5)
    ax.xaxis.set_major_formatter(hfmt)
    ax.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best',frameon=False)
    plt.xticks(rotation=30)
    plt.show()
    return

def rs_plot(fds,bi,marker=',',linestyle='solid'):
    """
    plota todos los 'adcval' de la estrutura
    'rs' contra el vector-tiempo contenido en fds
    """
    hfmt = dates.DateFormatter('%H:%M')
    fig = plt.figure(figsize=(11.69*0.75,8.27*0.75),num=None, dpi=100,facecolor='w',edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(fds,bi['adcval'],marker=marker,linestyle=linestyle)
    ax.xaxis.set_major_formatter(hfmt)
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_minor_locator(dates.MinuteLocator([10,20,30,40,50]))
    ax.xaxis.grid(True,'major',lw=2)
    ax.xaxis.grid(True,'minor',lw=0.5)
    plt.legend([b1,b2,b3,b4,b5,b6], ["beam 1","beam 2","beam 3","beam 4","beam 5" ,"beam 6"], loc='best',frameon=False)
    plt.xticks(rotation=30)
    plt.show()
    return

def bx_plot(taxdict,ybot,ytop):
    """
    Boxplots del contenido del dictionario (taxdict)
    devuelto por la funcion 'taxonomia'. Se tiene que
    introducir los limites del eje vertical: 'ybot', 'ytop'.
    """
    fig = plt.figure(figsize=(11.69,8.27), facecolor='w', edgecolor='k')
    ax = fig.add_subplot(191)
    ax.set_title('cold/hot')
    ax.boxplot(taxdict['colval'],0,'b+')#ax.scatter(coldfds,coldadc[:,0],color='blue',marker='.')
    ax.boxplot(taxdict['hotval'],0,'r+')#ax.scatter(hotfds,hotadc[:,0],color='red',marker='.')
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(192)
    ax.set_title('trk')
    ax.boxplot(taxdict['trkval'],0,'g+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(193)
    ax.set_title('map')
    ax.boxplot(taxdict['mapval'],0,'c+')#ax.scatter(taufds,tauadc[:,0],color='skyblue')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(194)
    ax.set_title('mpi')
    ax.boxplot(taxdict['mpival'],0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(195)
    ax.set_title('scn')
    ax.boxplot(taxdict['scnval'],0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(196)
    ax.set_title('sci')
    ax.boxplot(taxdict['scival'],0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(197)
    ax.set_title('tau')
    ax.boxplot(taxdict['tauval'],0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(198)
    ax.set_title('stl')
    ax.boxplot(taxdict['stlval'] ,0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    ax = fig.add_subplot(199)
    ax.set_title('unk')
    ax.boxplot(taxdict['unkval'],0,'y+')#ax.scatter(scanazfds,scanazadc[:,0],color='yellow')
    [label.set_visible(False) for label in ax.get_yticklabels()]
    plt.ylim(ybot,ytop)
    plt.grid()
    plt.show()
    return

def tx_plot(taxdict,str_key,boto=0,topo=7000,marker=',',linestyle='solid'):
    """
    Plot de los adcval de cualquiera de las variables segregadas
    segun el diccionario (taxdict) devuelto por la funcion 'taxonomia'.
    Se tiene que entrar con la clave (str_key) para escoger el modo, la
    misma que aparece en los titulos del grafico porducido por 'caja'.
    """
    strfds = str_key+'time'
    stradc = str_key+'val'
    hfmt = dates.DateFormatter('%H:%M:%S')
    fig = plt.figure(figsize=(11.69*0.75,8.27*0.75),num=None, dpi=100,facecolor='w',edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    [b1,b2,b3,b4,b5,b6] = ax.plot(taxdict[strfds],taxdict[stradc][:,0:6],marker=marker,linestyle=linestyle)
    plt.legend([b1,b2,b3,b4,b5,b6], ["1","2","3","4","5","6"], loc='best',frameon=False)
    ax.xaxis.set_major_formatter(hfmt)
    ax.xaxis.set_major_locator(dates.HourLocator())
    ax.xaxis.set_minor_locator(dates.MinuteLocator([10,20,30,40,50]))
    ax.xaxis.grid(True,'major',lw=2)
    ax.xaxis.grid(True,'minor',lw=0.5)
    plt.xticks(rotation=0)
    plt.ylim(bottom=boto,top=topo)
    plt.show()
    return

def st_plot(statdict):
    """
    Plot de la relacion lineal entre las estadisticas
    de cada canal de 212 GHz en relacion al canal 1 y
    del canal 5 (405 GHz) en relacion al canal 6
    """
    fig = plt.figure(figsize=(12,6),num=None,dpi=100,facecolor='w',edgecolor='k')
    ax1 = fig.add_subplot(1,2,1)
    ax1.plot(statdict['stat01'],statdict['stat01'],'k',ls='--',label = '1 vs 1')
    ax1.plot(statdict['stat01'],statdict['stat02'],'g+',label = '1 vs 2')
    ax1.plot(statdict['stat01'],statdict['stat03'],'r>',label='1 vs 3')
    ax1.plot(statdict['stat01'],statdict['stat04'],'cv',label='1 vs 4')
    plt.legend(loc='best',frameon=False)
    plt.xticks(rotation=30)
    #plt.grid()

    ax2 = fig.add_subplot(1,2,2)
    ax2.plot(statdict['stat06'],statdict['stat06'],'k',ls='--',label ='6 vs 6')
    ax2.plot(statdict['stat06'],statdict['stat05'],'y*',label='6 vs 5')
    plt.legend(loc='best',frameon=False)
    plt.xticks(rotation=30)
    #plt.grid()
    plt.show()
    return

def field_plot(pos,mr,fi,leftxlim=-1200,rightxlim=1200,bottomylim=-1200,topylim=1200):
    '''
    Plot del campo que contiene al sol y la progresion del centro
    de los feixes durante el periodo observado

    '''
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    circle212 = plt.Circle((0,0),mr,color='b',fill=False)
    circle405 = plt.Circle((0,0),mr,color='r',fill=False)
    ax.add_artist(circle212)
    ax.add_artist(circle405)
    ax.plot(pos.beams_ew[fi,:],pos.beams_ns[fi,:])
    ax.plot(pos.source_ew[fi],pos.source_ns[fi])
    ax.set_xlim(leftxlim,rightxlim)
    ax.set_ylim(bottomylim,topylim)
    plt.grid()
    return

def oldfield_plot(pos,mr,fi,leftxlim=-1200,rightxlim=1200,bottomylim=-1200,topylim=1200):
    '''
    Plot del campo que contiene al sol y la progresion del centro
    de los feixes durante el periodo observado
    '''
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    circle212 = plt.Circle((0,0),mr,color='b',fill=False)
    circle405 = plt.Circle((0,0),mr,color='r',fill=False)
    ax.add_artist(circle212)
    ax.add_artist(circle405)
    for j in range(0,6):
        ax.plot(2*pos.source_ew[fi]-pos.beams_ew[fi,j],2*pos.source_ns[fi]-pos.beams_ns[fi,j])
    ax.plot(pos.source_ew[fi],pos.source_ns[fi])
    ax.set_xlim(leftxlim,rightxlim)
    ax.set_ylim(bottomylim,topylim)
    plt.grid()
    return


def bgfit_plot(Xplusone,rshh,iniz,fina,dX,hc=False,left=0,right=1,bottom=10,top=3600):
    '''


    '''
    X = Xplusone-1
    obsadcX = rshh['adcval'][:,X][iniz:fina]# - np.min(rshh['adcval'][:,0][iniz:fina]).astype(float)
    fig = plt.figure(figsize=(11.69,8.27),num=None, dpi=100, facecolor='w', edgecolor='k')
    model_curve = obsadcX-dX ##- np.nanmin(obsadcX-dX)
    ax = fig.add_subplot(111)
    ax.plot(rshh['time']/3.6e7,model_curve,label='model map beam'+np.str(Xplusone),linewidth = 2)
    ax.plot(rshh['time']/3.6e7,obsadcX,label='obser map beam'+np.str(Xplusone),linewidth = 1.3)
    #ax.plot(rs.time/3.6e7,d1,label='diff mod-obs 1',linewidth = 0.62)
    #ax.plot(rs.time/3.6e7, predict_y, 'k-')
    #ax.plot(rs.time/3.6e7, pred_error, 'y-')
    plt.xlim(left=left,right=right)
    plt.ylim(bottom=bottom,top=top)
    ax.legend(loc='best')
    plt.grid()
    plt.show
    #print residual_std_error1
    if hc==True:
        savefig('/home/fer/Desktop/'+'sellsmoke.eps', papertype = 'a4', orientation = 'landscape', format = 'eps')
    return

def hat_plot(a,xoff,yoff,zi,w,h,myrstep=40,mycstep=40,myalpha=0.45,eleva=45,azimu=120,asksave=False):
    '''
    '''
    ssl = np.size(a[:,0])/2
    wlong = np.linspace(-ssl,ssl-1,w)
    hlong = np.linspace(-ssl,ssl-1,h)
    xfield,yfield=np.meshgrid(wlong,hlong)
    fig = plt.figure(figsize=(8.27*1.25,8.27*1.25),num=None,dpi=100,facecolor='w',edgecolor='k')
    ax = fig.add_subplot(1,1,1, projection='3d')

    ax.plot_surface(xfield[ssl-500:ssl+499,ssl-500:ssl+499],yfield[ssl-500:ssl+499,ssl-500:ssl+499],a[ssl-500:ssl+499,ssl-500:ssl+499],\
                rstride=myrstep, cstride=mycstep, alpha=myalpha,cmap=cm.hot)
    ax.plot_wireframe(xoff,yoff,zi,color='k',lw=1.)
    ax.set_xlabel('Azimuth offset (arcmin)')
    ax.set_ylabel('Elevation offset (arcmin)')
    ax.set_zlabel('Temperature (K)')
    #ax.set_xticklabels([30,20,10,0,-10,-20,-30])
    #ax.set_yticklabels([30,20,10,0,-10,-20,-30])
    axmin = 3.6*np.array(ax.get_xticks().tolist())/60.
    aymin = 3.6*np.array(ax.get_yticks().tolist())/60.
    ax.set_xticklabels(axmin)
    ax.set_yticklabels(aymin)
    #for angle in range(0, 360):
    #    ax.view_init(30, angle)
    #   plt.draw()
    ax.view_init(eleva,azimu)
    plt.grid()
    #plt.show()
    if asksave == True:
        plt.savefig('/home/fer/Desktop/'+'solar_map.eps', papertype = 'a4', orientation = 'portrait', format = 'eps')
    return
