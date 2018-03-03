from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d
from astropy import units as u
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.io import readsav
from scipy.interpolate import interp1d
from scipy import arange, array, exp
from scipy.interpolate import InterpolatedUnivariateSpline

def dered_extinc(wave,mag,color):
	wave = 1/wave # inverse microns
	
	inv_lam = [0.00,0.377,0.820,1.667,1.828,2.141,2.433,3.0704,3.846]
	extn_ratio = [0.00,0.265,0.829,2.688,3.055,3.806,4.315,6.625,6.591]

	f = interp1d(inv_lam, extn_ratio, kind='cubic')
	
	print 'org: ',mag,'f: ',f(1/wave),'new: ',mag+f(1/wave)*color

	return mag+f(1/wave)*color

def mag2mJy(mag,f_zpt):
	return f_zpt*10**(mag/(-2.5))*1000


def sed_catalog_search(target):

	#print target
	wave = np.array([])*u.nm
	phot = np.array([])*u.mJy
	phot_e = np.array([])*u.mJy
	notes = np.array([])

	#Thompson TD1 UV Catalog ----------

	cat = 'II/59B'
	v = Vizier(columns=['F2740','e_F2740','F2365','e_F2365','F1965','e_F1965','F1565','e_F1565'])
	res = v.query_object(target,catalog=[cat])

	if not not res:
		TD1_2740 = res[0]['F2740'][0]
		TD1_2740_e = res[0]['e_F2740'][0]

		TD1_2365 = res[0]['F2365'][0]
		TD1_2365_e = res[0]['e_F2365'][0]

		TD1_1965 = res[0]['F1965'][0]
		TD1_1965_e = res[0]['e_F1965'][0]

		TD1_1565 = res[0]['F1565'][0]
		TD1_1565_e = res[0]['e_F1565'][0]

		print "TD1"
		print TD1_2740,TD1_2365,TD1_1965,TD1_1565
		print TD1_2740_e,TD1_2365_e,TD1_1965_e,TD1_1565_e

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

	cat = 'I/259'
	v = Vizier(columns=['BTmag','e_BTmag','VTmag','e_VTmag'])
	res = v.query_object(target,catalog=[cat])

	if not not res:
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

		tychoB_e = np.abs(tychoB - tychoB_e)
		tychoV_e = np.abs(tychoV - tychoV_e)

		print "Tycho"
		print tychoV,tychoB
		print tychoB_e,tychoV_e

		wave = np.append(wave,np.array([4400,5500])*0.0001) #angstroms to microns
		phot = np.append(phot,[tychoB,tychoV])
		phot_e = np.append(phot_e,[tychoB_e,tychoV_e])
		notes = np.append(notes,['Tycho B','Tycho V'])

	#ASAA of ROSAT Catalog --------------------

	cat = 'J/AcA/62/67'
	v = Vizier(columns=['__Vmag_','__Imag_','sVmag','sImag'])
	res = v.query_object(target,catalog=[cat])

	if not not res:
		asasV = res[0]['__Vmag_'][0]
		asasI = res[0]['__Imag_'][0]

		asasV_e = res[0]['sVmag'][0]	
		asasI_e = res[0]['sImag'][0]

		asasV_e = mag2mJy(asasV_e+asasV,3781)
		asasI_e = mag2mJy(asasI_e+asasI,2635)

		asasV = mag2mJy(asasV,3781)	
		asasI = mag2mJy(asasI,2635)

		asasV_e = np.abs(asasV - asasV_e)
		asasI_e = np.abs(asasI - asasI_e)

		print "ROSAT"
		print asasV,asasI
		print asasV_e,asasI_e

		wave = np.append(wave,np.array([5500,9700])*0.0001) #angstroms to microns
		phot = np.append(phot,[asasV,asasI])
		phot_e = np.append(phot_e,[asasV_e,asasI_e])
		notes = np.append(notes,['ASAS V','ASAS I'])

	# 2MASS Catalog --------------------

	cat = 'II/246'
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

		print 'J-H',twomassJ-twomassH
		print 'H-K',twomassH-twomassK

		twomassJ_e = mag2mJy(twomassJ_e+twomassJ,1594)
		twomassH_e = mag2mJy(twomassH_e+twomassH,1024)
		twomassK_e = mag2mJy(twomassK_e+twomassK,666.7)

		twomassJ = mag2mJy(twomassJ,1594)
		twomassH = mag2mJy(twomassH,1024)
		twomassK = mag2mJy(twomassK,666.7)

		twomassJ_e = np.abs(twomassJ - twomassJ_e)
		twomassH_e = np.abs(twomassH - twomassH_e)
		twomassK_e = np.abs(twomassK - twomassK_e)

		print "2MASS"
		print twomassJ,twomassH,twomassK
		print twomassJ_e,twomassH_e,twomassK_e

		wave = np.append(wave,[1.235,1.662,2.159])
		phot = np.append(phot,[twomassJ,twomassH,twomassK])
		phot_e = np.append(phot_e,[twomassJ_e,twomassH_e,twomassK_e])
		notes = np.append(notes,['2MASS J','2MASS H','2MASS K'])


	#Spitzer Catalog --------------------

	cat = 'J/ApJS/211/25/catalog' # Chen 2014 Catalog; units microns and mJy
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

		print "Spitzer"
		print sptzr13irs,sptzr24mips,sptzr31irs,sptzr70mips
		print sptzr13irs_e,sptzr24mips_e,sptzr31irs_e,sptzr70mips_e

		wave = np.append(wave,[10.75,24,32,70])
		phot = np.append(phot,[sptzr13irs,sptzr24mips,sptzr31irs,sptzr70mips])
		phot_e = np.append(phot_e,[sptzr13irs_e,sptzr24mips_e,sptzr31irs_e,sptzr70mips_e])
		notes = np.append(notes,['Spitzer IRS 13','Spitzer MIPS 24','Spitzer IRS 31','Spitzer MIPS 70'])

	#WISE Catalog --------------------------

	cat ='II/328/allwise'
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

		print "WISE"
		print w1mag,w2mag,w3mag,w4mag
		print w1mag_e,w2mag_e,w3mag_e,w4mag_e

		wave = np.append(wave,[3.3526,4.6028,11.5608,22.0883])
		phot = np.append(phot,[w1mag,w2mag,w3mag,w4mag])
		phot_e  = np.append(phot_e,[w1mag_e,w2mag_e,w3mag_e,w4mag_e])
		notes = np.append(notes,['WISE 1','WISE 2','WISE 3','WISE 4'])

	#AKARI data
	cat='II/297'
	v = Vizier(columns=['S09','e_S09','S18','e_S18'])
	res = v.query_object(target,catalog=[cat],radius=Angle(3, "arcsec"))

	if not not res:
		#Jy to mJy
		S09mag = res[0]['S09'][0]*1000
		S18mag = res[0]['S18'][0]*1000

		S09mag_e = np.float(res[0]['e_S09'][0])*1000
		S18mag_e = np.float(res[0]['e_S18'][0])*1000

		if ((S09mag > 0) and (S09mag_e > 0)): 
			wave = np.append(wave,[8.61])
			phot = np.append(phot,[S09mag])
			phot_e = np.append(phot_e,[S09mag_e])
			notes = np.append(notes,['AKARI S9'])

		if ((S18mag > 0) and (S18mag_e > 0)): 
			wave = np.append(wave,[18.39])
			phot = np.append(phot,[S18mag])
			phot_e = np.append(phot_e,[S18mag_e])
			notes = np.append(notes,['AKARI L18'])

		print "AKARI"
		print S09mag,S18mag
		print S09mag_e,S18mag_e

	

	# END catalogs -------------------------------------

	return [wave, phot, phot_e, notes]

def pow_law(x,exp,C):
	return 10**(exp*np.log10(x)+C)

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
	    x=np.log(x)
	    y=np.log(y)
	    xs2=np.log(xs)
	    ys2=np.log(ys)
            return 10**(ys2[0]+(x-xs2[0])*(ys2[1]-ys2[0])/(xs2[1]-xs2[0]))
        elif x > xs[-1]:
	    x=np.log(x)
	    y=np.log(y)
	    xs2=np.log(xs)
	    ys2=np.log(ys)
            return 10**(ys2[-1]+(x-xs2[-1])*(ys2[-1]-ys2[-2])/(xs2[-1]-xs2[-2]))
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

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

	return [wave, phot, phot_e, notes]

def fit_sed(target,result,fit):

	wave, phot, phot_e, notes = list(result)

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
		f[id] = f[id]*(lam[id]/lam_o)**beta

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
		f1[id] = f1[id]*(lam[id]/lam_o)**beta
		f2[id] = f2[id]*(lam[id]/lam_o)**beta

		return f1+f2


	if fit.lower() == "single":
		#plot bb fit
		x1 = np.arange(0.2,4000,0.05)

		x0 = [5500,1]
		id = np.where(wave < 10)
		starpopt, pcov = curve_fit(sed_bb,wave[id],phot[id],sigma=phot_e[id],p0=x0)

		#starpopt = [5500,45]
		f1 = sed_bb(x1,*starpopt)

		#x0 = [300,1*10**6,150,1*10**7,210,1]
		x0 = [50,1*10**6,210.0,0.0]
		id = np.where(wave > 24)
		starbb = sed_bb(wave,*starpopt)

		dustpopt, pcov = curve_fit(mod_sed_bb,wave[id],phot[id]-starbb[id],sigma=phot_e[id],p0=x0)
		f2 = mod_sed_bb(x1,*dustpopt)
		diskbb = mod_sed_bb(wave,*dustpopt)

		#residual points
		resid = np.abs(phot-starbb-diskbb)
		ids = np.where(resid < 3*phot_e)
		plt.loglog(wave[ids],resid[ids],'yv',markersize=5)
		ids = np.where(resid > 3*phot_e)
		plt.loglog(wave[ids],resid[ids],'yo',markersize=5)

		ftot = f1+f2
		#plt.loglog(x_m,f*x_m**2,'b-')
		plt.loglog(x1,f1,'r--')
		plt.loglog(x1,f2,'g--')
		plt.loglog(x1,ftot,'k-')
		plt.loglog(x1,mod_sed_bb(x1,dustpopt[0],dustpopt[1],210,-1.0),'k--')
		plt.loglog(x1,mod_sed_bb(x1,dustpopt[0],dustpopt[1],210,-2.0),'k-.')

		plt.loglog(wave,phot,'co')
		plt.title(target)
		plt.ylim([10**(-1),10**4])
		plt.xlim([3*10**(-1),2*10**3])

		plt.xlabel('microns')
		plt.ylabel('mJy')

		print 'Stellar Teff: '+str(starpopt[0])
		print 'Dust Teff: '+str(dustpopt[0])

		print "Output File: "+str(target+'_'+fit+'.png')
		plt.savefig(target+'_'+fit+'.png')

	elif fit.lower() == "double":
	
		#plot bb fit
		x1 = np.arange(0.01,40000,0.05,dtype=np.float64)

		x0 = [5500,35]
		id = np.where(wave < 1)
		starpopt, pcov = curve_fit(sed_bb,wave[id],phot[id],sigma=phot_e[id],p0=x0)

		#starpopt = [5500,45]
		f1 = sed_bb(x1,*starpopt)

		x0 = [100,1*10**6,60,1*10**7,120,0]
		id = np.where(wave > 11)
		starbb = sed_bb(wave,*starpopt)

		dustpopt, dustpcov = curve_fit(dub_mod_sed_bb,wave[id],phot[id]-starbb[id],sigma=phot_e[id],p0=x0,maxfev=3000)
		f2 = dub_mod_sed_bb(x1,*dustpopt)

		diskbb = dub_mod_sed_bb(wave,*dustpopt)

		#residual points

		resid = np.abs(phot-starbb-diskbb)
		ids = np.where(resid < 3*phot_e)

		plt.loglog(wave[ids],resid[ids],'yv',markersize=5)
		ids = np.where(resid > 3*phot_e)
		plt.loglog(wave[ids],resid[ids],'yo',markersize=5)

		resid = np.abs(phot-starbb)
	
		plt.loglog(wave,resid,'yo',markersize=2)

		ftot = f1+f2

		plt.loglog(x1,f1,'r--')

		diskbb1 = mod_sed_bb(x1,dustpopt[2],dustpopt[3],dustpopt[4],dustpopt[5])
		diskbb2 = mod_sed_bb(x1,dustpopt[0],dustpopt[1],dustpopt[4],dustpopt[5])

		plt.loglog(x1,diskbb1,'g--')
		plt.loglog(x1,diskbb2,'g--')

		plt.loglog(x1,ftot,'k-',markersize=2)
		plt.xlabel(r'wavelength ($\mu$m)',fontsize="20")
		plt.ylabel('F (mJy)',fontsize="20")

		print 'Stellar Teff: '+str(starpopt[0])
		print 'Dust Teff: '+str(dustpopt[1])
		print 'Dust Teff: '+str(dustpopt[3])

		print "Output File: "+str(target+'_'+fit+'.png')
		plt.savefig(target+'_'+fit+'.png')


	
