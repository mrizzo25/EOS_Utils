#Monica Rizzo 2017
#
#this script interfaces with main integrator to 
#produce useful plots
#
#user inputs desired options on command-line
#Ex:python Make_Plots.py --eos-directory "HW_Glendenning5.10" --mass-radius

import numpy as np
import argparse
from scipy.interpolate import spline
import TOV
import matplotlib.pyplot as plt
import os

#scaling constants
msun=1.989*10**33
rho_nuc=2.3*10**14

#integration defaults
r_init=0.00001
r_step_init=0.0001
r_step_min=10**-15 #don't want to go lower than machine precision...
frac_error=0.0005
shooting_steps=200
shooting_tol=100


#command line arguments to generate various plots
parser=argparse.ArgumentParser()
parser.add_argument("--eos-directory",type=str,help="directory containing eos to use")
parser.add_argument("--mass-radius", default=None, action='store_true', help="generate mass-radius curve and related plots")
parser.add_argument("--save-table", action="store_true", help="save mass radius table if generating curve")
parser.add_argument("--stability-curves", default=None, action='store_true', help="plot a mass-radius stability plot from given directory")
parser.add_argument("--use-radius",type=float, help="use shooting method to find a given radius *not guarunteed to be physical*")
opts=parser.parse_args()

if opts.eos_directory:
   eos_dir=opts.eos_directory
else:
   print("Please specify EOS directory")

if opts.mass_radius:
   
   #scale central density
   scale=np.linspace(1.0,20.0,100)

   mass_array=np.zeros(len(scale))
   radius_array=np.zeros(len(scale))

   #generate plots of mass vs. radius for each starting central density
   i=0
   plt.figure(1)
   for s in scale:
      pm_array,r_array,rho_array=TOV.integrate_TOV(eos_dir,10.0,r_init,r_step_init,r_step_min,s*rho_nuc,frac_error,"pressure")
      plt.plot(r_array*10**-5,pm_array[:,1]/msun,'-o',label="rho_c="+str(s*rho_nuc))
      mass_array[i]=pm_array[-1,1]
      radius_array[i]=r_array[-1]
      print("Central density:   "+str(s*rho_nuc)+" g/cm^3    Final Mass: "+str(mass_array[i]/msun)+" M_sun    Final Radius: "+str(radius_array[i]*10**-5)+" km")
      i+=1


   if opts.save_table:  
       #save table and put in MR directory
       np.savetxt(str(eos_dir)+"_mr.dat", np.c_[scale*rho_nuc,radius_array*10**-5,mass_array/msun],header="rho_c   radius (km)  mass (m_sun)")
       cmd="mv "+str(eos_dir)+"_mr.dat MR/"
       os.system(cmd)

   #plot mass-radius curve
   plt.figure(2)
   plt.plot(radius_array*10**-5, mass_array/msun, '-bo')
   plt.xlabel("Mass (M_sun)")
   plt.ylabel("Radius (km)")
   plt.savefig(eos_dir+"_mass_rad.eps")

 
if opts.stability_curves:
   scale=np.array([10**(-9),10**(-8),10**(-6),10**(-4),0.1,0.2,0.3,0.4,0.5,0.85,0.9,1.0,1.5,2.0,2.5,3.0,4.0,5.0])

   mass_array=np.zeros(len(scale))
   radius_array=np.zeros(len(scale))

   #generate stability curve using each starting central density
   i=0
   plt.figure(5)
   for s in scale:
       pm_array,r_array,rho_array=TOV.integrate_TOV(eos_dir,10.0,r_init,r_step_init,r_step_min,s*rho_nuc,frac_error,"pressure")
       mass_array[i]=pm_array[-1,1]
       radius_array[i]=r_array[-1]
       i+=1

   plt.semilogx(radius_array*10**(-5),mass_array/msun,'-o')
   plt.xlabel("Radius (km)")
   plt.ylabel("Mass (M_sun)")
   plt.savefig(eos_dir+"_stability_regions.eps")


   plt.figure(6)
   plt.semilogx(scale*rho_nuc,mass_array/msun,'-o')
   plt.xlabel("Central Density (g/cm^3)")
   plt.ylabel("Mass (M_sun)")
   plt.savefig(eos_dir+"_density_mass_stability.eps")

if opts.use_radius:
   r_fin=opts.use_radius

   pm_array,r_array,rho_array=TOV.TOV_Shooting(eos_dir,r_fin,rho_nuc,shooting_steps,shooting_tol)

   plt.figure(7)
   plt.semilogx(r_array*10**-5,pm_array[:,1]/msun,'-o')
   plt.xlabel("Radius (km)")
   plt.ylabel("Mass (M_sun)")

   plt.savefig(eos_dir+"_mvr_shooting_radius_"+str(r_fin)+".eps")


