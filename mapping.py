#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit
import sys
sys.path.append("/home/fer/Documents/myPYTHONstuff/modules/")
from mapping import *
from other import *



def campo(tsky,tsun,solar_radius):
    '''
    Cria un campo do ceu de 6000x6000 pixels (6x6 degraus)
    e temperattura tsky em cujo centro cria um Sol artificial
    de radio solar_radius e temperatura tsun. 
    '''
    #dimensiones del campo (3.6 arcsec per pixel o 1 miligrado por pixel)
    w = np.int(6000) 
    h = w # the same for height
    #el radio del sol artificial escalado al campo
    #factor=1.0
    #innerradius = np.int(factor*mir/3.6)
    radius212 = np.int(solar_radius/3.6)
    #creando el sol artificial 212 y 405
    field = tsky*np.ones((w, h)).astype('float') # the background
    #the circle
    cx, cy = np.int(w/2), np.int(h/2)
    x_frame, y_frame = np.ogrid[-radius212:radius212, -radius212:radius212]
    index_mask = x_frame**2 + y_frame**2 <= radius212**2
    #index_central = np.where(x_frame**2 + y_frame**2 <= innerradius**2)
    field[cx-radius212:cx+radius212, cy-radius212:cy+radius212][index_mask] = tsun
    return field

def campolimbo(tsky,tsun,solar_radius):
    ''' 
    Cria un campo do ceu de 6000x6000 pixels (6x6 degraus)
    e temperattura tsky em cujo centro cria um Sol artificial
    de radio solar_radius e temperatura tsun. 
    '''
    #dimensiones del campo (3.6 arcsec per pixel o 1 miligrado por pixel)
    w = np.int(6000) 
    h = w # the same for height
    #el radio del sol artificial escalado al campo
    #factor=1.0
    #innerradius = np.int(factor*mir/3.6)
    radius212 = np.int(solar_radius/3.6)
    #creando el sol artificial 212 y 405
    field = tsky*np.ones((w, h)).astype('float') # the background
    #the circle
    cx, cy = np.int(w/2), np.int(h/2)
    x_frame, y_frame = np.ogrid[-radius212:radius212, -radius212:radius212]
    index_mask = x_frame**2 + y_frame**2 <= radius212**2
    #index_central = np.where(x_frame**2 + y_frame**2 <= innerradius**2)
    field[cx-radius212:cx+radius212, cy-radius212:cy+radius212][index_mask] = tsun

    lb0=1.0000*tsun
    ##the limb brightening
    hrad = np.int(1.0*radius212)
    irad = np.int(0.999*radius212)
    hcx, hcy = np.int(w/2), np.int(h/2)
    hy_frame, hx_frame = np.ogrid[-hrad: hrad, -hrad: hrad]
    hindex_mask = (hx_frame**2 + hy_frame**2 >= irad**2) & (hx_frame**2 + hy_frame**2 <= hrad**2)
    field[hcy-hrad:hcy+hrad, hcx-hrad:hcx+hrad][hindex_mask] = tsun+lb0

    lb1=0.20*tsun
    ##the limb brightening
    hrad = np.int(0.995*radius212)
    irad = np.int(0.85*radius212)
    hcx, hcy = np.int(w/2), np.int(h/2)
    hy_frame, hx_frame = np.ogrid[-hrad: hrad, -hrad: hrad]
    hindex_mask = (hx_frame**2 + hy_frame**2 >= irad**2) & (hx_frame**2 + hy_frame**2 <= hrad**2)
    field[hcy-hrad:hcy+hrad, hcx-hrad:hcx+hrad][hindex_mask] = tsun+lb1

    lb2=0.15*tsun
    ##the limb brightening
    hrad = np.int(0.85*radius212)
    irad = np.int(0.65*radius212)
    hcx, hcy = np.int(w/2), np.int(h/2)
    hy_frame, hx_frame = np.ogrid[-hrad: hrad, -hrad: hrad]
    hindex_mask = (hx_frame**2 + hy_frame**2 >= irad**2) & (hx_frame**2 + hy_frame**2 <= hrad**2)
    field[hcy-hrad:hcy+hrad, hcx-hrad:hcx+hrad][hindex_mask] = tsun+lb2

    lb3=0.10*tsun
    ##the limb brightening
    hrad = np.int(0.65*radius212)
    irad = np.int(0.35*radius212)
    hcx, hcy = np.int(w/2), np.int(h/2)
    hy_frame, hx_frame = np.ogrid[-hrad: hrad, -hrad: hrad]
    hindex_mask = (hx_frame**2 + hy_frame**2 >= irad**2) & (hx_frame**2 + hy_frame**2 <= hrad**2)
    field[hcy-hrad:hcy+hrad, hcx-hrad:hcx+hrad][hindex_mask] = tsun+lb3

    lb4=0.05*tsun
    ##the limb brightening
    hrad = np.int(0.35*radius212)
    irad = np.int(0.00*radius212)
    hcx, hcy = np.int(w/2), np.int(h/2)
    hy_frame, hx_frame = np.ogrid[-hrad: hrad, -hrad: hrad]
    hindex_mask = (hx_frame**2 + hy_frame**2 >= irad**2) & (hx_frame**2 + hy_frame**2 <= hrad**2)
    field[hcy-hrad:hcy+hrad, hcx-hrad:hcx+hrad][hindex_mask] = tsun+lb4

    ##the AR
    ar0 = 0.15*tsun
    arrad = np.int(np.ceil(150/3.6))
    arposx,arposy = -484/3.6, 203/3.6##-0.9*230,-0.5*230
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar0+tsun

    ##the AR
    ar1 = 0.15*tsun
    arrad = np.int(np.ceil(150./3.6))
    arposx,arposy = -254/3.6,165/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar1+tsun

    ar2 = 0.15*tsun
    arrad = np.int(np.ceil(150/3.6))
    arposx,arposy = 257/3.6, 99/3.6##-0.9*230,-0.5*230
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar2+tsun

    ##the AR
    ar3 = 0.15*tsun
    arrad = np.int(np.ceil(150./3.6))
    arposx,arposy = -475/3.6,-363/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar3+tsun

    ## post-flare region
    ##the AR
    ar4 = 0.15*tsun
    arrad = np.int(np.ceil(50./3.6))
    arposx,arposy = -757/3.6,-160/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar4+tsun


    ar5 = 0.15*tsun
    arrad = np.int(np.ceil(150/3.6))
    arposx,arposy = 845/3.6, 126/3.6##-0.9*230,-0.5*230
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar5+tsun

    ##the AR
    ar6 = 0.15*tsun
    arrad = np.int(np.ceil(150./3.6))
    arposx,arposy = 784/3.6,170/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar6+tsun

    ar7 = 0.15*tsun
    arrad = np.int(np.ceil(150/3.6))
    arposx,arposy = -392/3.6, -225/3.6##-0.9*230,-0.5*230
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar7+tsun

    ##the AR
    ar8 = 0.15*tsun
    arrad = np.int(np.ceil(150./3.6))
    arposx,arposy = 291/3.6,-146/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar8+tsun

    ## post-flare region
    ##the AR
    ar9 = 0.15*tsun
    arrad = np.int(np.ceil(50./3.6))
    arposx,arposy = 931/3.6,-146/3.6#-0.78*radius212,-0.63*radius212
    arcx, arcy = np.int(w/2)+arposx, np.int(h/2)+arposy
    ary_frame, arx_frame = np.ogrid[-arrad: arrad, -arrad: arrad]
    arindex_mask = arx_frame**2 + ary_frame**2 <= arrad**2 
    field[arcx-arrad:arcx+arrad, arcy-arrad:arcy+arrad][arindex_mask] = ar9+tsun
    return field
