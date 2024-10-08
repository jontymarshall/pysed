from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d
from astropy import units as u
from astropy import constants as c
from astropy.io import ascii,fits,votable
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.io import readsav
from scipy.interpolate import interp1d
from scipy import arange, array, exp
from scipy.interpolate import InterpolatedUnivariateSpline
import copy
from extinction import apply, remove, fitzpatrick99

def dered_extinc(wave,flux,Rv=3.1,Av=1.0):
    
    dered = remove(fitzpatrick99(wave.to('AA'), Av, Rv,unit='aa'), flux)
    
    return dered

def mag2mJy(mag,f_zpt):
    return f_zpt*10**(mag/(-2.5))*1000

def sed_vizier_votable(vot_file):
    
    votd = votable.parse(vot_file).get_first_table()
    wavs =(c.c / (1e9*votd.array['_sed_freq'].data/u.s)).to('micron')
    
    pick = np.where((votd.array['_sed_eflux'].mask == False)&(votd.array['_sed_eflux'].data != 0.)&(wavs < 2499*u.micron))
    wave = wavs[pick]
    phot = votd.array['_sed_flux'].data[pick]*u.Jy
    phot_e = votd.array['_sed_eflux'].data[pick]*u.Jy
    notes = votd.array['_sed_filter'].data[pick]
    
    fin = np.where(np.isfinite(phot_e))
    sed_data = [wave[fin], phot[fin], phot_e[fin], notes[fin]]
    
    return sed_data

def sed_catalog_search(target,catalogues=['II/59B','I/259','J/AcA/62/67','II/246','J/ApJS/211/25/catalog','II/328/allwise','II/297','VIII/106/hppsc070','VIII/106/hppsc100','VIII/106/hppsc160']):
    
    #list of Vizier catalogues which will return a query.
    master_catalogues = ['I/355/gaiadr3','II/349/ps1','II/293/glimpse','II/59B','I/259','J/AcA/62/67','II/246','J/ApJS/211/25/catalog','II/328/allwise','II/297','VIII/106/hppsc070','VIII/106/hppsc100','VIII/106/hppsc160']
    
    #print target
    wave = np.array([])*u.nm
    phot = np.array([])*u.mJy
    phot_e = np.array([])*u.mJy
    notes = np.array([])

    #Thompson TD1 UV Catalog ----------
    for cat in catalogues:
        
        if cat not in master_catalogues:
            print("Catalogue ",cat," not recognized; please check input or update master_catalogues.")
        
        else: 
            if cat == 'I/355/gaiadr3':
                v = Vizier(columns=['Gmag','e_Gmag','BPmag','e_BPmag','RPmag','e_RPmag'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    
                    fbp = res[0]['BPmag'][0]
                    e_fbp = res[0]['e_BPmag'][0]
                    
                    fg = res[0]['Gmag'][0]
                    e_fg = res[0]['e_Gmag'][0]
                    
                    frp = res[0]['RPmag'][0]
                    e_frp = res[0]['e_RPmag'][0]
                    
                    fbp_e = mag2mJy(fbp+e_fbp,3552)
                    fg_e  = mag2mJy(fg+e_fg,3296)
                    frp_e = mag2mJy(frp+e_frp,2555)
                    
                    fbp = mag2mJy(fbp,3552)
                    fg  = mag2mJy(fg,3296)
                    frp = mag2mJy(frp,2555)
                    
                    u_fbp = abs(fbp - fbp_e)
                    u_fg  = abs(fg - fg_e)
                    u_frp = abs(frp - frp_e)
                    
                    u_wave = ([504, 637, 762]*u.nm).to('micron')
                    
                    wave = np.append(wave,u_wave)
                    phot = np.append(phot,[fbp,fg,frp]*u.mJy)
                    phot_e = np.append(phot_e,[u_fbp,u_fg,u_frp]*u.mJy)
                    notes = np.append(notes,['Gaia BP','Gaia G','Gaia RP'])
            
            if cat == 'II/293/glimpse':
                v = Vizier(columns=['F(3.6)','e_F(3.6)','F(4.5)','e_F(4.5)','F(5.8)','e_F(5.8)','F(8.0)','e_F(8.0)'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                
                if not not res:
                    
                    f3p6 = res[0]['F_3.6_'][0]
                    e_f3p6 = res[0]['e_F_3.6_'][0]
                    
                    f4p5 = res[0]['F_4.5_'][0]
                    e_f4p5 = res[0]['e_F_4.5_'][0]
                    
                    f5p8 = res[0]['F_5.8_'][0]
                    e_f5p8 = res[0]['e_F_5.8_'][0]
                    
                    f8p0 = res[0]['F_8.0_'][0]
                    e_f8p0 = res[0]['e_F_8.0_'][0]
                    
                    u_wave = ([3600, 4500, 5800, 8000,]*u.nm).to('micron')
                    
                    wave = np.append(wave,u_wave)
                    phot = np.append(phot,[f3p6,f4p5,f5p8,f8p0]*u.mJy)
                    phot_e = np.append(phot_e,[e_f3p6,e_f4p5,e_f5p8,e_f8p0]*u.mJy)
                    notes = np.append(notes,['GLIMPSE 3p6','GLIMPSE 4p5','GLIMPSE 5p8','GLIMPSE 8p0'])
                    
            if cat == 'II/59B':
                v = Vizier(columns=['F2740','e_F2740','F2365','e_F2365','F1965','e_F1965','F1565','e_F1565'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    
                    TD1_2740 = res[0]['F2740'][0]
                    TD1_2740_e = res[0]['e_F2740'][0]
                    
                    TD1_2365 = res[0]['F2365'][0]
                    TD1_2365_e = res[0]['e_F2365'][0]
                    
                    TD1_1965 = res[0]['F1965'][0]
                    TD1_1965_e = res[0]['e_F1965'][0]
                    
                    TD1_1565 = res[0]['F1565'][0]
                    TD1_1565_e = res[0]['e_F1565'][0]
                    
                    # print("TD1")
                    # print(TD1_2740,TD1_2365,TD1_1965,TD1_1565)
                    # print(TD1_2740_e,TD1_2365_e,TD1_1965_e,TD1_1565_e)
                    
                    u_wave = ([2740, 2365, 1965, 1565,]*u.AA).to('micron')
                    u_phot = np.array([TD1_2740,TD1_2365,TD1_1965,TD1_1565])*1e-2 #centiWatts to Watts
                    u_phot_e = np.array([TD1_2740_e,TD1_2365_e,TD1_1965_e,TD1_1565_e])*1e-2 #centiWatts to Watts
                    
                    u_phot = u_phot * (u.Watt/u.meter**2/u.nm).to('mJy', equivalencies= u.spectral_density(u_wave) )
                    u_phot_e = u_phot_e * (u.Watt/u.meter**2/u.nm).to('mJy', equivalencies= u.spectral_density(u_wave) )
                    
                    wave = np.append(wave,u_wave)
                    phot = np.append(phot,u_phot)
                    phot_e = np.append(phot_e,u_phot_e)
                    notes = np.append(notes,['TD 1','TD 1','TD 1','TD 1'])
            
            #Tycho Catalog --------------------
            if cat == 'I/259':
                v = Vizier(columns=['BTmag','e_BTmag','VTmag','e_VTmag'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    tycho_wave = ([4400,5500,]*u.AA).to('micron')
                    
                    tychoB = res[0]['BTmag'][0]
                    tychoV = res[0]['VTmag'][0]
                    BmV = tychoB - tychoV
                    
                    tychoV = tychoV-0.090*(BmV)
                    tychoB = 0.850*(BmV)+tychoV
                    
                    color = BmV-0.438
                    
                    tychoB = dered_extinc(4400*0.0001,tychoB,color)
                    tychoV = dered_extinc(5500*0.0001,tychoV,color)
                    
                    tychoB_e = res[0]['e_BTmag'][0]
                    tychoV_e = res[0]['e_VTmag'][0]
                    
                    tychoB_e = mag2mJy(tychoB+tychoB_e,4130)
                    tychoV_e = mag2mJy(tychoV+tychoV_e,3781)
                    
                    tychoB = mag2mJy(tychoB,4130)
                    tychoV = mag2mJy(tychoV,3781)
                    
                    tychoB_e = abs(tychoB - tychoB_e)
                    tychoV_e = abs(tychoV - tychoV_e)
                    
                    # print("Tycho")
                    # print(tychoV,tychoB)
                    # print(tychoB_e,tychoV_e)
                    
                    wave = np.append(wave,tycho_wave) #angstroms to microns
                    phot = np.append(phot,[tychoB,tychoV]*u.mJy)
                    phot_e = np.append(phot_e,[tychoB_e,tychoV_e]*u.mJy)
                    notes = np.append(notes,['Tycho B','Tycho V'])
            
            #ASAA of ROSAT Catalog --------------------
            if cat == 'J/AcA/62/67':
                v = Vizier(columns=['__Vmag_','__Imag_','sVmag','sImag'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    asas_wave = ([5500,9700,]*u.AA).to('micron')
                    
                    asasV = res[0]['__Vmag_'][0]
                    asasI = res[0]['__Imag_'][0]
                    
                    asasV_e = res[0]['sVmag'][0]
                    asasI_e = res[0]['sImag'][0]
                    
                    asasV_e = mag2mJy(asasV_e+asasV,3781)
                    asasI_e = mag2mJy(asasI_e+asasI,2635)
                    
                    asasV = mag2mJy(asasV,3781)
                    asasI = mag2mJy(asasI,2635)
                    
                    asasV_e = abs(asasV - asasV_e)
                    asasI_e = abs(asasI - asasI_e)
                    
                    # print("ROSAT")
                    # print(asasV,asasI)
                    # print(asasV_e,asasI_e)
                    
                    wave = np.append(wave,asas_wave) #angstroms to microns
                    phot = np.append(phot,[asasV,asasI]*u.mJy)
                    phot_e = np.append(phot_e,[asasV_e,asasI_e]*u.mJy)
                    notes = np.append(notes,['ASAS V','ASAS I'])
            
            # 2MASS Catalog --------------------
            if cat == 'II/246':
                v = Vizier(columns=['Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    twomassJ = res[0]['Jmag'][0]
                    twomassH = res[0]['Hmag'][0]
                    twomassK = res[0]['Kmag'][0]
                    
                    twomassJ_e = res[0]['e_Jmag'][0]
                    twomassH_e = res[0]['e_Hmag'][0]
                    twomassK_e = res[0]['e_Kmag'][0]
                    
                    #twomassJ = dered_extinc(1.235,twomassJ,color)
                    #twomassH = dered_extinc(1.662,twomassH,color)
                    #twomassK = dered_extinc(2.159,twomassK,color)
                    
                    # print('J-H',twomassJ-twomassH)
                    # print('H-K',twomassH-twomassK)
                    
                    twomassJ_e = mag2mJy(twomassJ_e+twomassJ,1594)
                    twomassH_e = mag2mJy(twomassH_e+twomassH,1024)
                    twomassK_e = mag2mJy(twomassK_e+twomassK,666.7)
                    
                    twomassJ = mag2mJy(twomassJ,1594)
                    twomassH = mag2mJy(twomassH,1024)
                    twomassK = mag2mJy(twomassK,666.7)
                    
                    twomassJ_e = abs(twomassJ - twomassJ_e)
                    twomassH_e = abs(twomassH - twomassH_e)
                    twomassK_e = abs(twomassK - twomassK_e)
                    
                    # print("2MASS")
                    # print(twomassJ,twomassH,twomassK)
                    # print(twomassJ_e,twomassH_e,twomassK_e)
                    
                    wave = np.append(wave,[1.235,1.662,2.159]*u.micron)
                    phot = np.append(phot,[twomassJ,twomassH,twomassK]*u.mJy)
                    phot_e = np.append(phot_e,[twomassJ_e,twomassH_e,twomassK_e]*u.mJy)
                    notes = np.append(notes,['2MASS J','2MASS H','2MASS K'])
            
            #Spitzer Catalog --------------------
            if cat == 'J/ApJS/211/25/catalog': # Chen 2014 Catalog; units microns and mJy
                v = Vizier(columns=['F13','e_F13','F24','e_F24','F31','e_F31','F70','e_F70'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    sptzr13irs = res[0]['F13'][0] # mJy to Jy
                    sptzr24mips = res[0]['F24'][0]
                    sptzr31irs = res[0]['F31'][0]
                    sptzr70mips = res[0]['F70'][0]
                    
                    sptzr13irs_e = res[0]['e_F13'][0]
                    sptzr24mips_e = res[0]['e_F24'][0]
                    sptzr31irs_e = res[0]['e_F31'][0]
                    sptzr70mips_e = res[0]['e_F70'][0]
                    
                    # print("Spitzer")
                    # print(sptzr13irs,sptzr24mips,sptzr31irs,sptzr70mips)
                    # print(sptzr13irs_e,sptzr24mips_e,sptzr31irs_e,sptzr70mips_e)
                    
                    wave = np.append(wave,[10.75,24,32,70]*u.micron)
                    phot = np.append(phot,[sptzr13irs,sptzr24mips,sptzr31irs,sptzr70mips]*u.mJy)
                    phot_e = np.append(phot_e,[sptzr13irs_e,sptzr24mips_e,sptzr31irs_e,sptzr70mips_e]*u.mJy)
                    notes = np.append(notes,['Spitzer IRS 13','Spitzer MIPS 24','Spitzer IRS 31','Spitzer MIPS 70'])
            
            #WISE Catalog --------------------------
            if cat == 'II/328/allwise':
                v = Vizier(columns=['W1mag','e_W1mag','W2mag','e_W2mag','W3mag','e_W3mag','W4mag','e_W4mag'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    w1mag = res[0]['W1mag'][0]
                    w2mag = res[0]['W2mag'][0]
                    w3mag = res[0]['W3mag'][0]
                    w4mag = res[0]['W4mag'][0]
                    
                    w1mag = mag2mJy(w1mag,309.540)
                    w2mag = mag2mJy(w2mag,171.787)
                    w3mag = mag2mJy(w3mag,31.674)
                    w4mag = mag2mJy(w4mag,8.363)
                    
                    w1mag_e = res[0]['e_W1mag'][0]
                    w2mag_e = res[0]['e_W2mag'][0]
                    w3mag_e = res[0]['e_W3mag'][0]
                    w4mag_e = res[0]['e_W4mag'][0]
                    
                    # print("WISE")
                    # print(w1mag,w2mag,w3mag,w4mag)
                    # print(w1mag_e,w2mag_e,w3mag_e,w4mag_e)
                    
                    wave = np.append(wave,[3.3526,4.6028,11.5608,22.0883]*u.micron)
                    phot = np.append(phot,[w1mag,w2mag,w3mag,w4mag]*u.mJy)
                    phot_e  = np.append(phot_e,[w1mag_e,w2mag_e,w3mag_e,w4mag_e]*u.mJy)
                    notes = np.append(notes,['WISE 1','WISE 2','WISE 3','WISE 4'])
            
            #AKARI data --------------------------
            if cat == 'II/297':
                v = Vizier(columns=['S09','e_S09','S18','e_S18'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
            
                if not not res:
                    #Jy to mJy
                    S09mag = (res[0]['S09'][0]*u.Jy).to('mJy')
                    S18mag = (res[0]['S18'][0]*u.Jy).to('mJy')
                    
                    S09mag_e = (np.float(res[0]['e_S09'][0])*u.Jy).to('mJy')
                    S18mag_e = (np.float(res[0]['e_S18'][0])*u.Jy).to('mJy')
                    
                    if ((S09mag > 0) and (S09mag_e > 0)): 
                        wave = np.append(wave,[8.61]*u.micron)
                        phot = np.append(phot,[S09mag])
                        phot_e = np.append(phot_e,[S09mag_e])
                        notes = np.append(notes,['AKARI S9'])
                    
                    if ((S18mag > 0) and (S18mag_e > 0)): 
                        wave = np.append(wave,[18.39]*u.micron)
                        phot = np.append(phot,[S18mag])
                        phot_e = np.append(phot_e,[S18mag_e])
                        notes = np.append(notes,['AKARI L18'])
                
                # print("AKARI")
                # print(S09mag,S18mag)
                # print(S09mag_e,S18mag_e)
            
            #Herschel Point Source Catalogue data 
            if cat == 'VIII/106/hppsc070':
                v = Vizier(columns=['Flux','rms'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    herp070   = res[0]['Flux'][0]
                    herp070_e = res[0]['rms'][0]
                    
                    wave = np.append(wave,[70.0]*u.micron)
                    phot = np.append(phot,[herp070]*u.mJy)
                    phot_e  = np.append(phot_e,[herp070_e]*u.mJy)
                    notes = np.append(notes,['HPACS70'])
            
            #Herschel Point Source Catalogue data 
            if cat == 'VIII/106/hppsc100':
                v = Vizier(columns=['Flux','rms'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    herp100   = res[0]['Flux'][0]
                    herp100_e = res[0]['rms'][0]
                    
                    wave = np.append(wave,[100.0]*u.micron)
                    phot = np.append(phot,[herp100]*u.mJy)
                    phot_e  = np.append(phot_e,[herp100_e]*u.mJy)
                    notes = np.append(notes,['HPACS100'])
            
            #Herschel Point Source Catalogue data 
            if cat == 'VIII/106/hppsc160':
                v = Vizier(columns=['Flux','rms'])
                res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))
                
                if not not res:
                    herp160   = res[0]['Flux'][0]
                    herp160_e = res[0]['rms'][0]
                    
                    wave = np.append(wave,[160.0]*u.micron)
                    phot = np.append(phot,[herp160]*u.mJy)
                    phot_e  = np.append(phot_e,[herp160_e]*u.mJy)
                    notes = np.append(notes,['HPACS160'])
            
    # END catalogs -------------------------------------
    
    # remove nan values
    fin = np.where(np.isfinite(phot_e))
    sed_data = [wave[fin], phot[fin], phot_e[fin], notes[fin]]
    print(np.transpose(sed_data))

    return sed_data

def pow_law(x,exp,C):
    return 10**(-2*np.log10(x)+C)

def extrap1d(xo,yo,xs,ys):
    
    xwav = copy.copy(xs)
    xflx = copy.copy(ys)
    
    if np.min(xo) < np.min(xs):
        
        extra_wave = np.logspace(np.log10(0.9*np.min(xo)),np.log10(xs[0]),num=1000,base=10.0,endpoint=True)
        
        slope_real = (ys[1] - ys[0]) / (xs[1] - xs[0])
        extra_real = ys[0] + (slope_real * (extra_wave - xs[0]))
        
        xwav = np.append(extra_wave,xwav[1:])
        xflx = np.append(extra_real,xflx[1:])
    
    if np.max(xo) > np.max(xs):
        
        extra_wave = np.logspace(np.log10(xs[-1]),np.log10(1.1*np.max(xo)),num=1000,base=10.0,endpoint=True)
        
        slope = (ys[-2] - ys[-1]) / (xs[-2] - xs[-1])
        extra_real = ys[-1] + (slope * (extra_wave - xs[-1]))
        
        xwav = np.append(xwav,extra_wave[1:])
        xflx = np.append(xflx,extra_real[1:])
    
    xwav = xwav*u.AA
    xflx = xflx*u.mJy
    
    return xwav, xflx

def read_spitzer_irsa(file):
    archive = zipfile.ZipFile(file,'r')
    files = zipfile.ZipFile.namelist(archive)
    tables = [s for s in files if 'tune.tbl' in s]
    
    data = [0,0,0,0,0]

    for t in tables:
        tabletoread = zipfile.ZipFile.read(archive,t)
        lines = tabletoread.split('\n')
        for l in lines:
            if ('char' in l) or ('|' in l) or ('\\' in l):
                pass
            else:
                dataline = filter( None,l.split('  ') )
                if (np.size(dataline) > 3):
                    data = np.vstack((data,dataline))
    
    return (data[1:,0:]).astype(np.float)

def add_photometry(result,phot_arr):
    
    wave, phot, phot_e, notes = list(result)
    
    wave = np.append(wave,phot_arr[0])
    phot = np.append(phot,phot_arr[1])
    phot_e = np.append(phot_e,phot_arr[2])
    notes = np.append(notes,phot_arr[3])
    
    noinf = np.isfinite(wave)
    
    wave = wave[noinf]
    phot = phot[noinf]
    phot_e = phot_e[noinf]
    notes = notes[noinf]
    
    return [wave, phot, phot_e, notes]

def fit_sed(target,result,nbb='single',star_type='bb',star_models={},distance=100,extinction=False):
    
    wave, phot, phot_e, notes = list(result)
    
    try:
        wave = wave.to('micron').value
    except:
        pass
    
    try:
        phot = phot.value
    except:
        pass
    
    try:
        phot_e = phot_e.value
    except:
        pass
    
    if extinction:
        phot = dered_extinc(wave*u.micron,phot,Rv=3.1,Av=distance/1000.) #Av ~ 1 mag/kpc in optical
    
    # print("Photometry for ",target)
    # print(wave,phot,phot_e)
    
    def sed_bb(lam,T,area):
        #lam in microns
        #T in Kelvin
    
        k=1.3806488*10**(-16) #boltzmann erg/K
        h=6.62606957*10**(-27) #plank erg*s
        c=2.99792458*10**14 #speed of light micron/s
    
        f=((2*h*c**2)/lam**5)/(np.exp((h*c)/(lam*k*T))-1) # F_nu in mJy
        f = f*lam**2
        
        return f*area
    
    def mod_sed_bb(lam,T,area,lam_o,beta):
        #lam in microns
        #T in Kelvin
        
        k=1.3806488*10**(-16) #boltzmann erg/K
        h=6.62606957*10**(-27) #plank erg*s
        c=2.99792458*10**14 #speed of light micron/s
        
        f=((2*h*c**2)/lam**5)*(1/(np.exp(h*c/(lam*k*T))-1)) # F_nu in mJy
        f = f*lam**2
        
        id = np.where(lam > lam_o)
        f[id] = f[id]*(lam_o/lam[id])**beta
        
        return f*area
    
    def dub_mod_sed_bb(lam,T1,area1,T2,area2,lam_o,beta):
        #lam in microns
        #T in Kelvin
        #lam_o=210.0
        #beta=0.8
        
        k=1.3806488*10**(-16) #boltzmann erg/K
        h=6.62606957*10**(-27) #plank erg*s
        c=2.99792458*10**14 #speed of light micron/s
        
        f1=((2*h*c**2)/lam**5)*(1/(np.exp(h*c/(lam*k*T1))-1)) # F_nu in mJy
        f1=f1*lam**2*area1
        
        f2=((2*h*c**2)/lam**5)*(1/(np.exp(h*c/(lam*k*T2))-1)) # F_nu in mJy
        f2=f2*lam**2*area2
        
        id = np.where(lam > lam_o)
        f1[id] = f1[id]*(lam_o/lam[id])**beta
        f2[id] = f2[id]*(lam_o/lam[id])**beta
        
        return f1+f2
    
    #--------------------
    
    x1 = np.arange(0.1,2500,0.05)
    
    fig, ax = plt.subplots(1,1,figsize=(8, 6),dpi=300)
    
    if star_type.lower() == "bb":
        #plot bb fit
        
        x0 = [5500,1]
        id = np.where(wave < 10)
        starpopt, pcov = curve_fit(sed_bb,wave[id],phot[id],sigma=phot_e[id],p0=x0)
        
        #starpopt = [5500,45]
        starbb = sed_bb(wave,*starpopt)
        fstar = sed_bb(x1,*starpopt)
        
        ax.loglog(x1,fstar,'r-')
        
        print('Stellar Teff: '+str(starpopt[0]))
        
    elif star_type.lower() == "stellar":
        #print(wave)
        id = np.where(wave < 10)
        #print(id)
        #print(wave[id],phot[id],phot_e[id])
        x_m,f_m,temp_str,logg_str = chisqr_stellar_models(wave[id],phot[id],phot_e[id],star_models,pin_wave = 1.5*u.micron)
        
        print('Stellar Teff  : '+temp_str)
        print('Stellar log(g): '+logg_str)
        
        starbb_func = interp1d(x_m,f_m)
        starbb = starbb_func(wave)
        
        try:
            starbb = starbb.value
        except: 
            pass
        
        fstar = starbb_func(x1)
        
        ax.loglog(x_m,f_m,'r-')
        
    if nbb.lower() == "single":

        #x0 = [300,1*10**6,150,1*10**7,210,1]
        x0 = [50,1e6,210.0,-2.0]
        id = np.where(wave >= 10)
        
        #print(id)
        #print(wave[id],phot[id],phot_e[id],starbb[id])
        #print(wave,phot,starbb,phot_e)
        
        #print("phot: ",phot)
        
        dustpopt, pcov = curve_fit(mod_sed_bb,wave[id],phot[id]-starbb[id],sigma=phot_e[id],p0=x0,bounds=([10.,1e-30,50.,-4.],[1000.,1e30,500.,4.]))
        
        f2 = mod_sed_bb(x1,*dustpopt)
        diskbb = mod_sed_bb(wave,*dustpopt)
        
        #residual points
        resid = np.abs(phot-starbb-diskbb)
        ids = np.where(resid < 3*phot_e)
        ax.loglog(wave[ids],resid[ids],'yv',markersize=5)
        ids = np.where(resid > 3*phot_e)
        ax.loglog(wave[ids],resid[ids],'yo',markersize=5)
        
        ftot = fstar+f2
        #ax.loglog(x_m,f*x_m**2,'b-')
        
        ax.loglog(x1,f2,'g--')
        ax.loglog(x1,ftot,'k-')
        #ax.loglog(x1,mod_sed_bb(x1,dustpopt[0],dustpopt[1],210,-1.0),'k--')
        #ax.loglog(x1,mod_sed_bb(x1,dustpopt[0],dustpopt[1],210,-2.0),'k-.')
        
        frac_lum = get_lum(x1,diskbb,distance)/get_lum(x1,ftot,distance)
        
        print('Dust Teff: '+str(dustpopt[0]))
        print('Modified BB Wavelength: '+"{:.2f}".format(dustpopt[2]))
        print('Modified BB Slope: '+"{:.2f}".format(dustpopt[3]))
        print('Log Fractional Luminosity: '+"{:.2f}".format(np.log10(frac_lum)))
        
        print("Output File: "+str(target+'_'+nbb+'.png'))
        
    elif nbb.lower() == "double":
        
        x0 = [200,1E6,50,1E7,210,0]
        id = np.where(wave > 10)
        
        bounds = ((5,0,5,0,120,-5),(2000,np.inf,2000,np.inf,500,5))
        
        dustpopt, dustpcov = curve_fit(dub_mod_sed_bb,wave[id],phot[id]-starbb[id],sigma=phot_e[id],p0=x0,bounds=bounds)
        f2 = dub_mod_sed_bb(x1,*dustpopt)
        
        diskbb = dub_mod_sed_bb(wave,*dustpopt)
        
        #residual points
        
        resid = np.abs(phot-starbb-diskbb)
        ids = np.where(resid < 3*phot_e)
        
        ax.loglog(wave[ids],resid[ids],'yv',markersize=5)
        ids = np.where(resid > 3*phot_e)
        ax.loglog(wave[ids],resid[ids],'yo',markersize=5)
        
        resid = np.abs(phot-starbb)
        
        ax.loglog(wave,resid,'yo',markersize=2)
        
        ftot = fstar+f2
        
        diskbb1 = mod_sed_bb(x1,dustpopt[2],dustpopt[3],dustpopt[4],dustpopt[5])
        diskbb2 = mod_sed_bb(x1,dustpopt[0],dustpopt[1],dustpopt[4],dustpopt[5])
        
        ax.loglog(x1,diskbb1,'g--')
        ax.loglog(x1,diskbb2,'g--')
        
        ax.loglog(x1,ftot,'k-',markersize=2)
        
        frac_lum = (get_lum(x1,diskbb1,distance)+get_lum(x1,diskbb2,distance))/get_lum(x1,ftot,distance)
        
        #print 'Stellar Teff: '+str(starpopt[0])
        print('Warm Dust Teff: '+"{:.2f}".format(dustpopt[0]))
        print('Cold Dust Teff: '+"{:.2f}".format(dustpopt[2]))
        print('Modified BB Wavelength: '+"{:.2f}".format(dustpopt[4]))
        print('Modified BB Slope: '+"{:.2f}".format(dustpopt[5]))
        print('Log Fractional Luminosity: '+"{:.2f}".format(np.log10(frac_lum)))
        #print dustpopt
        
        print("Output File: "+str(target+'_'+nbb+'.png'))
    
    elif nbb.lower() == "none":
        pass
    
    ax.loglog(wave,phot,'co')
    ax.set_xlabel(r'wavelength ($\mu$m)',fontsize="20")
    ax.set_ylabel('F (mJy)',fontsize="20")
    ax.set_title(target,fontsize="20")
    ax.set_ylim([np.min(phot)/10, np.max(phot)*10])
    ax.set_xlim([np.min(wave)/10, np.max(wave)*10])
    
    plt.savefig(target+'_'+nbb+'.png')

    return ax

def compile_stellar_models(folder,fixed_logg=False,logg=4.5):
    from astropy.io import fits
    import glob
    
    star_models = {}
    if fixed_logg == True:
        files = glob.glob(folder+'lte*'+str(logg)+'*.fits')
    else:
        files = glob.glob(folder+'lte*.fits')
    
    if (np.size(files) < 1):
        raise ValueError('No stellar models have been found. Go download them.')

    for f in files:
        #print(f)
        
        hdulist = fits.open(f)
        flx = hdulist[0].data
        x = hdulist[0].header['CRVAL1']+np.arange(hdulist[0].header['NAXIS1'])*hdulist[0].header['CDELT1']
        x_m = x.astype(np.float)*1e-4 #Ang to microns
        flx = np.sum(flx,axis=0)*((x_m*1e4)**2)/2.998E14 # convert to mJy
        hdulist.close()

        f2 = np.convolve(flx, np.ones((15,))/15)
        sz1 = np.shape(flx)[0]
        sz2 = np.shape(f2)[0]
        bound = int(np.abs(sz1-sz2))
        f2 = f2[bound:-bound-1]
        x_m = x_m[bound//2:-(bound//2)-1]
        #print np.size(x_m),np.size(f2)

        rjtw = np.logspace(np.log10(np.max(x_m)),np.log10(2500),100)

        rtpopt,rtpcov = curve_fit(pow_law,x_m[-2000:],f2[-2000:],p0=[2.00,9])

        rjtf = pow_law(rjtw,*rtpopt)

        f2 = np.append(f2,rjtf)
        x_m = np.append(x_m,rjtw)
        
        #plt.loglog(x_m,f2,'b-')
        #plt.show
        
        if '+' in f:
            model_temp = (str(f).replace(folder,'').split('+')[0]).replace('lte','')
            model_logg = (str(f).replace(folder,'').split('-')[0].split('+')[1])
        else:
            model_temp = (str(f).replace(folder,'').split('-')[0]).replace('lte','')
            model_logg = (str(f).replace(folder,'').split('-')[1]).replace('.PHOENIX','')
        
        temp = np.float(model_temp)
        logg = np.float(model_logg)
        
        star_models[str(temp),str(logg)] = np.vstack((x_m,f2))
        
    return star_models

def chisqr_stellar_models(star_wave,star_phot,star_phot_e,star_models,pin_wave = 1.5*u.micron):

    model_flx = np.arange(np.size(star_wave))
    temps = np.array([0])
    loggs = np.array([0])
    
    for star in star_models.items():
        
        temp,logg = star[0]
        x_m,f2 = star[1]
        
        func = interp1d(x_m,f2,kind='nearest')
        
        model_flx = np.vstack((model_flx,[func(star_wave)]))
        temps = np.vstack((temps,np.float(temp)))
        loggs = np.vstack((loggs,np.float(logg)))
    
    model_flx = model_flx[1:,0:]
    
    temps = temps[1:,0:]
    loggs = loggs[1:,0:]
    
    #print(star_wave,pin_wave)
    
    #print(np.min(abs(star_wave.to('micron').value - pin_wave.value)))
    
    pin_id = np.where(np.min(abs(star_wave - pin_wave.value)) == abs(star_wave - pin_wave.value))
    
    if len(pin_id[0]) != 1:
        pin_id = ([np.array(pin_id[0][0])],)
        
    
    #print(star_phot[pin_id])
    #print(model_flx[0:,pin_id])
    
    diff_flx = model_flx[0:,pin_id]/star_phot[pin_id]
    
    #print diff_flx[0:,0]
    
    diff_flx_arr = np.repeat(diff_flx[0:,0],np.shape(model_flx)[1],axis=1)
    
    star_phot_arr = np.repeat([star_phot],np.shape(model_flx)[0],axis=0)
    star_phot_e_arr = np.repeat([star_phot_e/star_phot**2],np.shape(model_flx)[0],axis=0)
    
    # print("Model flux ",model_flx)
    # print("Diff flux ",diff_flx.value)
    # print("Star flux ",star_phot_arr)
    
    #print(model_flx.shape,model_flx)
    #print(diff_flx_arr.shape,diff_flx_arr)
    #print(star_phot_arr.shape,star_phot_arr)
    
    #print((model_flx/diff_flx_arr)-star_phot_arr)
    
    chisqr = np.sum(((model_flx/diff_flx_arr-star_phot_arr)**2)*star_phot_e_arr,axis=1)
    
    min_id = np.where(np.min(chisqr) == chisqr)[0][0]
    temp_str = str(temps[min_id][0])
    logg_str = str(loggs[min_id][0])
    model_diff = diff_flx[min_id][0][0]
    
    x_m, f_m = star_models[temp_str,logg_str]
    
    return x_m,f_m/model_diff,temp_str,logg_str

def get_lum(x1,f1,distance):
    import numpy as np

    lightspeed = 2.99792458*10**14 #micron/s
    fq = lightspeed/x1**2
    f1b = (f1/10**29)*fq #(F_nu to F_lam)
    Jysum = np.trapz(x1,f1b)
    return (Jysum*4*np.pi*distance**2)/(3.826*10**26)

