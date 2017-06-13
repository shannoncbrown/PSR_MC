#This script grids over different parameters (D, h3 to start with) for PSR J1640+2224
import os
import numpy as np 
import subprocess #Might be more useful for running bash commands?
import astropy.units as u 
import astropy.units.astrophys as ua 
import astropy.constants as con 

#Define some early parameters
name_of_par_file = 'something.par' #The par file shouldn't contain the parameters we're gridding over, since we'll need to write these values to the temp file eventually
name_of_tim_file = 'something.tim'

#Define some new units for us to use
ls = ua.lyr * (1./31557600.00)
micros = u.s * 1.e-6
r_units = ua.solMass * micros

e = 0.00079727109617760987333 #from update.par file
c_con = con.c #m/s
pm_ra = 2.083752 * u.marcsec/u.yr #This is mult by cos_dec, mas/year
err_pm_ra =  0.064657 * u.marcsec/u.yr #mas/year
pm_dec = -11.041313 * u.marcsec/u.yr #mas/year
err_pm_dec = 0.097153 * u.marcsec/u.yr #mas/year
pm_tot = np.sqrt(pm_ra**2. + pm_dec**2.) #mas/year
theta_mu = np.arccos(pm_dec/pm_tot) #Don't need to convert from mas to rad because units cancel
Pb = 175.46066194 *u.d #days, from Lohmer
x = 55.329716380086932295 * ls#light-seconds, from update.par
err_x = 0.00000008056985585768 * ls
T_sun = 4.925490947 * micros #In microseconds
f = 0.0059074304 * ua.solMass #Mass function

#Define functions to convert between parameters
def s_func(cosi):
	return np.sqrt(-cosi**2. + 1.)

def r_func(h3, cosi):
	return h3 * ((1. + cosi)/(1. - cosi))**1.5

def x_dot_func(pm_tot, x, cosi, theta_mu, big_omega):
	return -x * pm_tot * (cosi/np.sin(abs(np.arccos(cosi)))) * np.sin(theta_mu - big_omega)

def d_func(px):
	return (1./px) * ua.kpc*u.marcsec

def Pb_dot_func(pm_ra, pm_dec, d, Pb):
	return (((pm_ra.to(u.rad/u.s))**2. + (pm_dec.to(u.rad/u.s) )**2.)*d.to(u.m)/c_con) * Pb.to(u.s)

def omega_dot_GR_func(M, e, Pb):
	return 3. * ((T_sun * M)**(2./3.)/(1. - e**2.)) * (Pb/(2.*np.pi))**(-5./3.)

def omega_dot_k_func(pm_tot, cosi, theta_mu, big_omega):
	return (pm_tot/np.sin(abs(np.arccos(cosi)))) * np.cos(theta_mu - big_omega)

def M_func(f, s, r):
	return np.sqrt((r/T_sun)**3. * s**3./f)

#Create a grid/cube of the parameters we want to do a MC for

px = np.arange(0.702548 - 0.060684,0.702548 +0.060684, 0.001) * u.marcsec #Taking the range of parallaxes within the error reported by VLBI measurements
h3 = np.arange(0.01, 10., 0.01) *r_units #Not sure what appropriate values would be for this
cosi = np.arange(0.01, 1, 0.01) #Dimensionless
big_omega = np.arange(0., 180., 5) *u.deg #Not sure what appropriate values for this are either, or if we'll just keep it constant

#Start a series of for loops
for p in px:
	for h in h3:
		for c in cosi:
			for o in big_omega:
				#Convert these values into values tempo2 will recognize using the functions defined about
				s = s_func(c)
			
				r = r_func(h, c)
				
				x_dot = x_dot_func(pm_tot, x, c, theta_mu, o) 
				print x_dot
				x_dot = x_dot.to(u.lyr*u.rad/u.s) 
				print x_dot
				d = d_func(p)
				Pb_dot = Pb_dot_func(pm_ra, pm_dec, d, Pb)
				print Pb_dot
				M = M_func(f, s, r) #calculate from m2 and i
				omega_dot_GR = omega_dot_GR_func(M.value, e, Pb)
				omega_dot_GR = omega_dot_GR.to(u.yr**(2./3.)/u.yr**(5./3.)) * 206264806.2471 * u.mas #Last factor is to convert to mas
				
				omega_dot_k = omega_dot_k_func(pm_tot, c, theta_mu, o)
				
				omega_dot = omega_dot_GR + omega_dot_k
				omega_dot = omega_dot.to(u.deg/u.yr) #To put it in the units needed by tempo2
				
				print omega_dot

				#Create a new parameter file with the MC parameters in it
				os.system('cp %s temp.par' %name_of_par_file)

				#Basic form for putting calculated parameters in .par file: os.system('echo \n "NAME OF VARIABLE IN TEMPO2		" %s "1" >> temp.par' %value_of_variable)
				os.system('echo "SINI		" %s "1" >> temp.par' %s)
				os.system('echo "M2		" %s "1" >> temp.par' %(r/T_sun.value))
				os.system('echo "PX		" %s "1" >> temp.par' %p.value)
				os.system('echo "XDOT		" %s "1" >> temp.par' %(x_dot.value*31557600.00)) #Need to convert back to lt-s/s
				os.system('echo "PBDOT		" %s "1" >> temp.par' %Pb_dot.value)
				os.system('echo "OMDOT		" %s "1" >> temp.par' %omega_dot.value)

				#Run tempo2 with this new par file, record the chi_squared
				#To grab a value from tempo2 James used the bash command $(/home/jmckee/tempo2runtime/bin/tempo2 -f temp.par $tim_file | grep "$fit_parameter" | awk '{print $5 " " $6}') 
				#How can we translate this into python?
				#Replace this with a bash script that does this and records the chi squared value, since python doesn't seem to like setting variables to shell output
				#chis = os.system('tempo2 -f temp.par %s | grep "TITLE OF CHI SQUARED IN TEMPO2 OUTPUT"' % name_of_tim_file) #I don't think we need the awk command, grep should be able to grab the value?


				os.system('bash mc_bash.sh %s %d %d %d' %(name_of_tim_file, p.value, h.value, c))

				quit()

