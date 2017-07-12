#Monica Rizzo, 2017
#Tolmann-Oppenheimer Volkoff (TOV) Equation Solver
#Needs: Tabulated equations of state (EOS), in same directory
#EOS files read and interpolated by EOS_Tables.py
#
#Solves TOV equation using pre-computed presssure-density. 
#Assumes simplified form where baryon number density is not used.
#
#Can be implemented as both a shooting method problem, and a simple
#forward integration problem
#
#Ex: TOV_Shooting('HW_Glendenning5.8',10.0,rho_nuc,200,10)
#takes an initial central density guess and integrates out to a certain radius;
#finds the central density which produces a desired radius
#NOTE: Some radii cannot be achieved for a given equation of state, so it 
# is useful to look at results of the forward integration method first
# to determine physical regimes
#Ex:integrate_TOV("BPS_Glendenning5.8",10.0,0.01,0.01,10**-16,rho_nuc,0.001,"pressure")
#integrates the TOV equation until the pressure reaches zero, for a given starting 
#central density (assumed to be around nuclear density)


import numpy as np
import EOSfromTable as eos
import matplotlib.pyplot as plt


#rk4 step for TOV
def TOV_rk4(pm,r,dr,rho,consts,p_array,rho_array,line_upper,line_lower):
    k1=TOV(r,pm,rho)
    

    rhok2=rho_of_p(pm[0]+0.5*k1[0]*dr,p_array,rho_array,consts,line_upper,line_lower)
    k2=TOV(r+0.5*dr,pm+0.5*k1*dr,rhok2)
    
    rhok3=rho_of_p(pm[0]+0.5*k2[0]*dr,p_array,rho_array,consts,line_upper,line_lower)
    k3=TOV(r+0.5*dr,pm+0.5*k2*dr,rhok3)
    
    rhok4=rho_of_p(pm[0]+k3[0]*dr,p_array,rho_array,consts,line_upper,line_lower)
    k4=TOV(r+dr,pm+k3*dr,rhok4)

    return pm+(1/6.0)*(k1+2*k2+2*k3+k4)*dr

#fundamental constants in cgs
G=6.67259*10**-8
sol=2.99792458*10**10
msun=1.989*10**33
rho_nuc=2.3*10**14



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
def TOV(r,pm,rho):
     P=pm[0]
     m=pm[1]
     dpdr= -G/r**2*(rho+P/sol**2)*(m+(4.*np.pi*r**3*P)/sol**2)*(1-(2.*G*m)/(r*sol**2))**(-1)
     dmdr= 4.*np.pi*r**2*rho
     return np.array([dpdr,dmdr])

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


#input radius in *CENTIMETERS*
def integrate_TOV(eos_name,r_fin,r_init,r_step_init,r_step_min,rhoc_guess,desired_fractional_error,termination):
   r_fin=r_fin*10**5
   e0=desired_fractional_error

   #import p and rho from table
   P_table,rho_table=eos.p_rho_arrays(eos_name)
  

   #compute constants for p(rho) and rho(p) functions
   p_const,p_line_upper,p_line_lower = eos.interp_eos_p_of_rho(eos_name)
   rho_const,rho_line_upper,rho_line_lower = eos.interp_eos_rho_of_p(eos_name)
   
   #initial conditions
   r = r_init
   rho=rhoc_guess
   p=p_of_rho(rho,P_table,rho_table,p_const,p_line_upper,p_line_lower)
   
   m=0.001
   pm_vec=np.array([p,m])  
   pm_array=np.array([pm_vec])
   r_array=np.array([r_init])
   rho_array=np.array([rhoc_guess]) 
   
   dr=r_step_init

   #print r_fin
   #for shooting method, terminate when a desired radius is reached
   if termination =="radius":
    while r<=r_fin:
     #two temporary vectors
     pm_vec1=pm_vec
     pm_vec2=pm_vec
 
     #twice a single step
     pm_vec1=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
   
     #one step, two times
     pm_vec2=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
     pm_vec2=TOV_rk4(pm_vec2,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)

     #calculate fractional error
     e_f=abs(pm_vec1-pm_vec2)/abs(pm_vec+TOV(r,pm_vec,rho)*dr)
     
     e_f=max(e_f[0],e_f[1])

     #keep step and increase
     if e_f<e0:
        if r+2*dr<=r_fin:
           r_array=np.append(r_array,r)
           pm_array=np.vstack((pm_array,pm_vec2-(pm_vec1-pm_vec2)/15.))
           rho_array=np.append(rho_array,rho_of_p((pm_vec2[0]-(pm_vec1[0]-pm_vec2[0])/15.),P_table,rho_table,rho_const,rho_line_upper,rho_line_lower))
        r=r+2*dr
        pm_vec=pm_vec2-(pm_vec1-pm_vec2)/15.
        rho=rho_of_p(pm_vec[0],P_table,rho_table,rho_const,rho_line_upper,rho_line_lower)
        #increase step
        if e_f==0.0:
           dr=15*dr
        elif 0.95*abs(e0/e_f)**(0.2)<5. and e_f!=0.0:
           dr=0.95*dr*abs(e0/e_f)**(0.2)
        else:
           dr=5*dr
        if r+2*dr>r_fin and r!=r_fin:
           dr=(r_fin-r)/2.
     if e_f>e0:
        if 0.95*dr*abs(e0/e_f)**(0.25)<r_step_min:
           print("minimum step size exceeded")
           return pm_array,r_array,rho_array
        else:
           dr=0.95*dr*abs(e0/e_f)**(0.25)
        if r+2*dr>r_fin and r!=r_fin:
           dr=(r_fin-r)/2.

   #for forward integration, terminate when the pressure drops to zero
   elif termination =="pressure":
    while pm_vec[0]>=10**(-15):
     #two temporary vectors
     pm_vec1=pm_vec
     pm_vec2=pm_vec
 
     #twice a single step
     pm_vec1=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
   
     #one step, two times
     pm_vec2=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
     pm_vec2=TOV_rk4(pm_vec2,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)

     #calculate fractional error
     e_f=abs(pm_vec1-pm_vec2)/abs(pm_vec+TOV(r,pm_vec,rho)*dr)
     
     e_f=max(e_f[0],e_f[1])

     #keep step and increase
     if e_f<e0:
        #if r+2*dr<=r_fin:
        r_array=np.append(r_array,r)
        pm_array=np.vstack((pm_array,pm_vec2-(pm_vec1-pm_vec2)/15.))
        rho_array=np.append(rho_array,rho_of_p((pm_vec2[0]-(pm_vec1[0]-pm_vec2[0])/15.),P_table,rho_table,rho_const,rho_line_upper,rho_line_lower))
        r=r+2*dr
        pm_vec=pm_vec2-(pm_vec1-pm_vec2)/15.
        rho=rho_of_p(pm_vec[0],P_table,rho_table,rho_const,rho_line_upper,rho_line_lower)
        #increase step
        if e_f==0.0:
           dr=15*dr
        elif 0.95*abs(e0/e_f)**(0.2)<5. and e_f!=0.0:
           dr=0.95*dr*abs(e0/e_f)**(0.2)
        else:
           dr=5*dr        
   
        boundary_overstepped=False
        try:
           pm_vec_tmp=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
           pm_vec_tmp1=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
           pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
        except:
           boundary_overstepped=True
       
        #adjust step at the boundary to not overshoot
        if boundary_overstepped==True:
           while boundary_overstepped==True or pm_vec_tmp[0]<0. or pm_vec_tmp1[0]<0.:
              dr=0.5*dr
              try:
                 pm_vec_tmp=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
                 pm_vec_tmp1=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
                 pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
              except:
                 boundary_overstepped=True
              else:
                 boundary_overstepped=False
              if (pm_vec_tmp1[0]>0. and pm_vec_tmp1[0]<1*10**-16):
                  break
              if (pm_vec_tmp[0]>0. and pm_vec_tmp[0]<1*10**-16):
                  break




     if e_f>e0:
        if 0.95*dr*abs(e0/e_f)**(0.25)<r_step_min:
           print("minimum step size exceeded")
           return pm_array,r_array,rho_array
        else:
           dr=0.95*dr*abs(e0/e_f)**(0.25)
   
        boundary_overstepped=False
        try:
           pm_vec_tmp=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
           pm_vec_tmp1=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
           pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)

        except:
           boundary_overstepped=True

        #adjust step at the boundary to not overshoot
        if boundary_overstepped==True:
           while boundary_overstepped==True or pm_vec_tmp[0]<0. or pm_vec_tmp1[0]<0.:
              dr=0.5*dr
              try:
                 pm_vec_tmp=TOV_rk4(pm_vec,r,2*dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
                 pm_vec_tmp1=TOV_rk4(pm_vec,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
                 pm_vec_tmp1=TOV_rk4(pm_vec_tmp1,r,dr,rho,rho_const,P_table,rho_table,rho_line_upper,rho_line_lower)
              except:
                 boundary_overstepped=True
              else:
                 boundary_overstepped=False
              if (pm_vec_tmp1[0]>0. and pm_vec_tmp1[0]<1*10**-16):
                  break
              if (pm_vec_tmp[0]>0. and pm_vec_tmp[0]<1*10**-16):
                  break



        #if p+2*dr>r_fin and r!=r_fin:
        #   dr=(r_fin-r)/2.
       
   return pm_array,r_array,rho_array


#pms,rs,rhos=integrate_TOV("HarrisonWheeler",0.1,1.0,0.01,0.001,rho_nuc*10**-6,0.001)
#plt.plot(rs*10**-5,(pms[:,1])/msun)
#plt.show()

#shooting method to solve TOV equation for a desired radius
def TOV_Shooting(eos_name,r_final,rho_c,n0,tol):
    P_table,rho_table=eos.p_rho_arrays(eos_name)
    a=rho_c*0.01
    if a<min(rho_table):
      a=0.2*min(rho_table)
    b=rho_c*10
    if b>max(rho_table):
       b=2*max(rho_table)
    top=b
    print b
    i=1
    while i<=n0:
       p=a+(b-a)/2.0
       p=p**10
       print p
       pm_array,r_array,rho_array=integrate_TOV(eos_name,r_final,0.01,0.01,10**-16,p,0.5,"radius")
       val=pm_array[-1,0]
       print("step: "+str(i)+"  mass: "+str(pm_array[-1,1]/msun)+"  radius:  "+str(r_array[-1]*10**-5)) 
       print (np.log(b)-np.log(a))/2.0
       if val==0 or (np.log(b)-np.log(a))/2.0<tol:
          if p==top:
                 print("Hitting top of table") 
          print("rho_c Value: "+str(p))
          return pm_array,r_array,rho_array
       i+=1
       pm_array,r_array,rho_array=integrate_TOV(eos_name,r_final,0.01,0.01,10**-16,a,0.5,"radius")
       val_a=pm_array[-1,0]
       if val_a*val>0:
          a=p
       else:
          b=p
    print("failed")
    return pm_array,r_array,rho_array


"""
testing block
radii=np.array([8.5,9.0,9.5,10.0,10.1,10.2,10.3,10.5,10.5,11.1])
mass_array=np.zeros(len(radii))
radius_array=np.zeros(len(radii))
i=0
f,ax=plt.subplots(2,2)
plt.figure(1)
for r in radii:
  pms,rs,rhos=TOV_Shooting('HW_Glendenning5.8',r,rho_nuc*10**-2,200,10)
 plt.plot(rs*10**-5,pms[:,1]/msun)
 plt.xlabel("Radius (km)")
 plt.ylabel("Mass (M_sun)")
 ax[0, 1].semilogy(rs*10**-5,pms[:,0],'o')
 ax[0, 1].set_xlabel("Radius (km)")
 ax[0, 1].set_ylabel("Pressure (dyne/cm^2)")
 ax[1, 0].semilogy(rs*10**-5,rhos,'o')
 ax[1, 0].set_xlabel("Radius (km)")
 ax[1, 0].set_ylabel("Density (g/cm^3)")
 ax[1, 1].loglog(rhos,pms[:,0],'o')
 ax[1, 1].set_xlabel("Density (g/cm^3)")
 ax[1, 1].set_ylabel("Pressure (dyne/cm^2)")

  mass_array[i]=pms[-1,1]
  radius_array[i]=rs[-1]
  i+=1


plt.savefig("HW_Glendenning5.8_plots.png")
plt.figure(2)
plt.plot(radius_array*10**-5,mass_array/msun,'-o')
plt.xlabel("Radius (km)")
plt.ylabel("Mass (M_sun)")
plt.savefig("HW_Glendenning5.8_mass_rad.png")

pms,rs,rhos=integrate_TOV("HW_Glendenning5.8",10.0,0.01,0.01,10**-16,rho_nuc,0.001)

testing interpolation
P_table,rho_table=eos.p_rho_arrays("HW_Glendenning5.8")
p_const,p_line_upper,p_line_lower = eos.interp_eos_p_of_rho("HW_Glendenning5.8")
rho_const,rho_line_upper,rho_line_lower = eos.interp_eos_rho_of_p("HW_Glendenning5.8")

rho_array=np.zeros(len(np.arange(5*10**30,1*10**33,10**30)))
p_array=np.arange(5*10**30,1*10**33,10**30)
i=0
for p in np.arange(5*10**30,1*10**33,10**30):
   rho_array[i]=rho_of_p(p_array[i],P_table,rho_table,rho_const,rho_line_upper,rho_line_lower)
   i+=1 

plt.figure(2)
plt.plot(p_array,rho_array,'o')
plt.plot(P_table[30:40],rho_table[30:40])
plt.show()

plt.figure(1)
plt.plot(rs*10**-5,(pms[:,1])/msun,'o')
plt.xlabel("Radius (km)")
plt.ylabel("Mass (M_sun)")
plt.figure(2)
plt.loglog(rs*10**-5,rhos)
plt.ylabel("Density (g/cm^3)")
plt.xlabel("Radius (km)")
plt.figure(3)
plt.loglog(rs*10**-5,pms[:,0])
plt.ylabel("Pressure (dyne/cm^2)")
plt.xlabel("Radius (km)")
plt.show()"""
