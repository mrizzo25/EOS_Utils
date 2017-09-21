#Monica Rizzo, 2017
#Tolmann-Oppenheimer Volkoff (TOV) Equation Solver
#Needs: Tabulated equations of state (EOS)
#EOS files read and interpolated by EOSfromTable.py
#
#Solves TOV equation using pre-computed pressure-density tables. 

import numpy as np
import EOSfromTable as eos
import matplotlib.pyplot as plt

#fundamental constants in cgs
G=6.67259*10**-8
sol=2.99792458*10**10
msun=1.989*10**33
rho_nuc=2.3*10**14


#rk4 step for TOV
def TOV_rk4(pm,r,dr,rho,nb,rho_consts,nb_consts,p_array,rho_array,nb_array,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower):
    k1=TOV(r,pm,rho,nb)

    rhok2=rho_of_p(pm[0]+0.5*k1[0]*dr,p_array,rho_array,rho_consts,rho_line_upper,rho_line_lower)
    nbk2=nb_of_p(pm[0]+0.5*k1[0]*dr,p_array,nb_array,nb_consts,nb_line_upper,nb_line_lower)
    k2=TOV(r+0.5*dr,pm+0.5*k1*dr,rhok2,nbk2)
    
    rhok3=rho_of_p(pm[0]+0.5*k2[0]*dr,p_array,rho_array,rho_consts,rho_line_upper,rho_line_lower)
    nbk3=nb_of_p(pm[0]+0.5*k2[0]*dr,p_array,nb_array,nb_consts,nb_line_upper,nb_line_lower)
    k3=TOV(r+0.5*dr,pm+0.5*k2*dr,rhok3,nbk3)
    
    rhok4=rho_of_p(pm[0]+k3[0]*dr,p_array,rho_array,rho_consts,rho_line_upper,rho_line_lower)
    nbk4=nb_of_p(pm[0]+k3[0]*dr,p_array,nb_array,nb_consts,nb_line_upper,nb_line_lower)
    k4=TOV(r+dr,pm+k3*dr,rhok4,nbk4)

    return pm+(1/6.0)*(k1+2*k2+2*k3+k4)*dr


def handle_boundary(pm,rad,dr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower):
  """
  Function to adjust step size near star boundary to prevent overstepping
  """
  delr=dr
  boundary_overstepped=False
  try:
         pm_vec_tmp=TOV_rk4(pm,rad,2*delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
         pm_vec_tmp1=TOV_rk4(pm,rad,delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
         pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,rad,delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)

  except:
         boundary_overstepped=True

  #adjust step at the boundary to not overshoot 
  if boundary_overstepped==True:
      while boundary_overstepped==True or pm_vec_tmp[0]<0. or pm_vec_tmp1[0]<0.:
         delr=0.5*delr
         try:
              pm_vec_tmp=TOV_rk4(pm,rad,2*delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
              pm_vec_tmp1=TOV_rk4(pm,rad,delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
              pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,rad,delr,rho,nb,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
              if (pm_vec_tmp1[0]>0. and pm_vec_tmp1[0]<1*10**-16):
                  return delr
              if (pm_vec_tmp[0]>0. and pm_vec_tmp[0]<1*10**-16):
                  return delr
         except:
           boundary_overstepped=True
         else:
           boundary_overstepped=False
  return delr



#integrate mass at a given radius assiming an equation of 
#state with density profile rho
def int_mass(r,rho,n):
    h = r/n
    s = rho + r**2*rho
    for i in np.arange(1,n-1,1):
       r=i*h
       if i%2 == 0:
         s+=2.0*r**2*rho
       else:
         s+=4.0*r**2*rho
    s = s*h/3.0*np.pi*4.0
    return s 



#TOV diff eq, returns dp/dr and dm/dr as vector
def TOV(r,pm,rho,bnd):
     P=pm[0]
     m=pm[1]
     dpdr= -G/r**2*(rho+P/sol**2)*(m+(4.*np.pi*r**3*P)/sol**2)*(1-(2.*G*m)/(r*sol**2))**(-1)
     dnbdr= 4.*np.pi*r**2*bnd/(1. - 2*G*m/(sol**2*r))**(1./2.)
     dmdr= 4.*np.pi*r**2*rho
     return np.array([dpdr,dmdr,dnbdr])

# find pressure at a given density using interpolated function
def p_of_rho(rho,p_array,rho_array,consts,line_upper,line_lower):
    rho=np.log10(rho)
    p_array=np.log10(p_array)
    rho_array=np.log10(rho_array)
    for i in np.arange(0,len(rho_array)-1,1):
       if rho_array[i]<rho and rho_array[i+1]>rho:
          return 10.**(p_array[i]+consts[i,0]*(rho-rho_array[i])+consts[i,1]*(rho-rho_array[i])**2+consts[i,2]*(rho-rho_array[i])**3)
       elif rho==rho_array[i]:
          return 10.**p_array[i] 
       elif rho==rho_array[i+1]:
          return 10**p_array[i+1] 
    if rho>rho_array[-1]:
       return 10.**(line_upper[0]*rho+line_upper[1])
    elif rho<rho_array[0]:
       return 10.**(line_lower[0]*rho+line_lower[1])
 

#find density at a given pressure using interpolated function
def rho_of_p(p,p_array,rho_array,consts,line_upper,line_lower):
    p=np.log10(p)
    p_array=np.log10(p_array)
    rho_array=np.log10(rho_array)
    for i in np.arange(0,len(p_array)-1,1):
       if p_array[i]<p and p_array[i+1]>p:
          return 10.**(rho_array[i]+consts[i,0]*(p-p_array[i])+consts[i,1]*(p-p_array[i])**2+consts[i,2]*(p-p_array[i])**3)
       elif p==p_array[i]:
          return 10.**rho_array[i]
       elif p==p_array[i+1]:
          return 10.**rho_array[i+1]
    if p>p_array[-1]:
       return 10.**(line_upper[0]*p+line_upper[1])
    elif p<p_array[0]:
       return 10.**(line_lower[0]*p+line_lower[1])


def nb_of_p(p,p_array,nb_array,consts,line_upper,line_lower):
    p=np.log10(p)
    p_array=np.log10(p_array)
    nb_array=np.log10(nb_array)
    for i in np.arange(0,len(p_array)-1,1):
       if p_array[i]<p and p_array[i+1]>p:
          return 10.**(nb_array[i]+consts[i,0]*(p-p_array[i])+consts[i,1]*(p-p_array[i])**2+consts[i,2]*(p-p_array[i])**3)
       elif p==p_array[i]:
          return 10.**nb_array[i]
       elif p==p_array[i+1]:
          return 10.**nb_array[i+1]
    if p>p_array[-1]:
       return 10.**(line_upper[0]*p+line_upper[1])
    elif p<p_array[0]:
       return 10.**(line_lower[0]*p+line_lower[1])


#input radius in *CENTIMETERS*
def integrate_TOV(eos_name,r_fin,r_init,r_step_init,r_step_min,rhoc_guess,desired_fractional_error,store_all=False):
   r_fin=r_fin*10**5
   e0=desired_fractional_error

   #import p and rho from table
   nb_table,P_table,rho_table=eos.p_rho_arrays(eos_name)
  

   #compute constants for p(rho) and rho(p) functions
   p_const,p_line_upper,p_line_lower = eos.interp_eos_p_of_rho(eos_name)
   rho_const,rho_line_upper,rho_line_lower = eos.interp_eos_rho_of_p(eos_name)
   nb_const,nb_line_upper,nb_line_lower = eos.interp_eos_nb_of_p(eos_name)
   
   #initial conditions
   r = r_init
   rho=rhoc_guess
   p=p_of_rho(rho,P_table,rho_table,p_const,p_line_upper,p_line_lower)
   bnd=nb_of_p(p,P_table,nb_table,nb_const,nb_line_upper,nb_line_lower)
   n=0.
   
   m=0.
   pm_vec=np.array([p,m,n])
   if store_all==True: 
     pm_array=np.array([pm_vec])
     r_array=np.array([r_init])
     rho_array=np.array([rhoc_guess]) 
   
   dr=r_step_init

   while pm_vec[0]>=10**(-15):
     #two temporary vectors
     pm_vec1=pm_vec
     pm_vec2=pm_vec
 
     #twice a single step
     pm_vec1=TOV_rk4(pm_vec,r,2*dr,rho,bnd,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
   
     #one step, two times
     pm_vec2=TOV_rk4(pm_vec,r,dr,rho,bnd,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)
     pm_vec2=TOV_rk4(pm_vec2,r,dr,rho,bnd,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)

     #calculate fractional error
     e_f=abs(pm_vec1-pm_vec2)/abs(pm_vec+TOV(r,pm_vec,rho,bnd)*dr)
     
     e_f=max(e_f[0],e_f[1])

     #keep step and increase
     if e_f<e0:
        #if r+2*dr<=r_fin:
        if store_all==True:
           r_array=np.append(r_array,r)
           pm_array=np.vstack((pm_array,pm_vec2-(pm_vec1-pm_vec2)/15.))
           rho_array=np.append(rho_array,rho_of_p((pm_vec2[0]-(pm_vec1[0]-pm_vec2[0])/15.),P_table,rho_table,rho_const,rho_line_upper,rho_line_lower))
        r=r+2*dr
        pm_vec=pm_vec2-(pm_vec1-pm_vec2)/15.
        rho=rho_of_p(pm_vec[0],P_table,rho_table,rho_const,rho_line_upper,rho_line_lower)
        bnd=nb_of_p(pm_vec[0],P_table,nb_table,nb_const,nb_line_upper,nb_line_lower)
        #increase step
        if e_f==0.0:
           dr=15*dr
        elif 0.95*abs(e0/e_f)**(0.2)<5. and e_f!=0.0:
           dr=0.95*dr*abs(e0/e_f)**(0.2)
        else:
           dr=5*dr        
   
        dr=handle_boundary(pm_vec,r,dr,rho,bnd,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower)


     if e_f>e0:
        if 0.95*dr*abs(e0/e_f)**(0.25)<r_step_min:
           if store_all==True:
               return pm_array,r_array,rho_array
           else:
               return pm_vec,r,rho
        else:
           dr=0.95*dr*abs(e0/e_f)**(0.25)
   
        dr=handle_boundary(pm_vec,r,dr,rho,bnd,rho_const,nb_const,P_table,rho_table,nb_table,rho_line_upper,rho_line_lower,nb_line_upper,nb_line_lower) 
      

        #if p+2*dr>r_fin and r!=r_fin:
        #   dr=(r_fin-r)/2.
   if store_all==True:
       return pm_array,r_array,rho_array
   else:
       return pm_vec,r,rho

