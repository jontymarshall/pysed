#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:52:26 2024

@author: jonty
"""

import pysed

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import units as u

phoenix_models = pysed.compile_stellar_models("/Users/jonty/mydata/phoenix/",fixed_logg=True,logg=4.50)

direc = "/Users/jonty/mydata/disks/hd105/"
#target_data = ascii.read(direc+"test_target_data.txt")

#targets   = target_data["name"].data
#distances = 1000./target_data["parallax"].data

targets   = ["HD 105"]
distances = [38.830] #pc

for i in range(len(targets)):
    
    target = targets[i]
    distance = distances[i]
    
    data = pysed.sed_catalog_search(target)
    
    if target == "HD 105":
        data = pysed.add_photometry(data,[160 *u.micron,112. *u.mJy,9. *u.mJy,"HPACS160"])
        data = pysed.add_photometry(data,[250 *u.micron,57. *u.mJy,10. *u.mJy,"HSPIRE250"])
        data = pysed.add_photometry(data,[350 *u.micron,38. *u.mJy,15. *u.mJy,"HSPIRE350"])
        data = pysed.add_photometry(data,[1270 *u.micron,2. *u.mJy,0.4 *u.mJy,"ALMA B6"])
    
    ax = pysed.fit_sed(target,data,nbb="double",star_type="stellar",star_models=phoenix_models,distance=distance)
    ax.plot()
    plt.show()
    plt.close()