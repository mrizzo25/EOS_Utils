#! /usr/bin/python

#Monica Rizzo 2017
#
#this script interfaces with main integrator to 
#produce useful plots and generate MR tables
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
baryon_mass=1.66*10**(-24)

#integration defaults
r_init=0.00001
r_step_init=0.0001
r_step_min=10**-17
frac_error=0.0005

#command line arguments to generate various plots
parser=argparse.ArgumentParser()
parser.add_argument("--eos-name",type=str,help="directory containing eos to use")
parser.add_argument("--mass-radius", default=None, action='store_true', help="generate mass-radius curve and related plots")
#parser.add_argument("--save-table", action="store_true", help="save mass radius table if generating curve")
#parser.add_argument("--stability-curves", default=None, action='store_true', help="plot a mass-radius stability plot from given directory")
parser.add_argument("--max-radius",type=float, default=18, help="define minimum radius to include in mass-radius curves")
parser.add_argument("--min-radius",default=None, help="if specified, radius to construct MR curve out to")
parser.add_argument("--save-plot", action='store_true', help="saves a plot of the mass radius curve generated")
opts=parser.parse_args()

if opts.eos_name:
   eos_dir=opts.eos_name
else:
   print("Please specify EOS directory")


fac_min=-0.2
if opts.mass_radius:

   r_fin=0
   while r_fin > 10+opts.max_radius or r_fin < opts.max_radius:
      print "Trying:"
      print "rho_c: ",(10**fac_min)*rho_nuc
      pm,r_fin,rho_fin=TOV.integrate_TOV(eos_dir,10.0,r_init,r_step_init,r_step_min,(10**fac_min)*rho_nuc,0.01)
      r_fin=r_fin*10**-5
      print "R: ",r_fin
      if r_fin<opts.max_radius:
         fac_min-=0.05
      elif r_fin>10+opts.max_radius:
         fac_min+=0.01

   fac_max=1.6
   if opts.min_radius!=None:
      r_fin=20.
      while r_fin>opts.min_radius or r_fin<opts.min_radius-2.:
          print "Trying:"
          print "rho_c: ",(10**fac_max)*rho_nuc
          pm,r_fin,rho_fin=TOV.integrate_TOV(eos_dir,10.0,r_init,r_step_init,r_step_min,(10**fac_max)*rho_nuc,0.01)
          r_fin=r_fin*10**-5
          print "R: ",r_fin
          if r_fin>opts.min_radius:
             fac_max+=0.05
          if r_fin<opts.min_radius-2.:
             fac_max-=0.01


   #scale central density
   scale=np.logspace(fac_min,fac_max)

   mass_array=np.zeros(len(scale))
   radius_array=np.zeros(len(scale))
   baryon_mass_array=np.zeros(len(scale))


   #generate mass vs. radius for each starting central density
   i=0
   for s in scale:
      pm_array,r,rho=TOV.integrate_TOV(eos_dir,10.0,r_init,r_step_init,r_step_min,s*rho_nuc,frac_error)
      mass_array[i]=pm_array[1]
      radius_array[i]=r
      baryon_mass_array[i]=pm_array[2]
      print("Central density:   "+str(s*rho_nuc)+" g/cm^3    Final Mass: "+str(mass_array[i]/msun)+" M_sun     Baryonic Mass:  "+str(baryon_mass_array[i]*baryon_mass/msun)+" M_sun    Final Radius: "+str(radius_array[i]*10**-5)+" km   Final Pressure: "+str(pm_array[0]))
      i+=1
   
   #check to make sure MR directory exists
   if "MR" not in os.listdir("./"):
      os.mkdir("MR")

   #save table and put in MR directory
   np.savetxt(str(eos_dir)+"_mr.dat", np.c_[scale*rho_nuc,radius_array*10**-5,mass_array/msun,baryon_mass_array*baryon_mass/msun],header="rho_c   radius (km)  mass (m_sun)  mb (m_sun)")
   cmd="mv "+str(eos_dir)+"_mr.dat MR/"
   os.system(cmd)

   if opts.save_plot:
      #plot mass-radius curve
      plt.figure(2)
      plt.plot(mass_array/msun,radius_array*10**-5,'-bo')
      plt.xlabel("Mass (M_sun)")
      plt.ylabel("Radius (km)")
      plt.savefig(eos_dir+"_mass_rad.eps")

 

