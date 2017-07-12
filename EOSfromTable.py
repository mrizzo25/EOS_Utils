import numpy as np 
import MonotonicSpline as ms

sol=2.99792458*10**10

#equations of state to choose from
eos = {'HW_Glendenning5.8':'/home/monica/Documents/School_Stuff/Ad_Comp/final_project/HW_Glendenning5.8/eos.dat','alf1':'/home/monica/Documents/REU_Materials/EOS_Tables/alf1.dat','sly':'/home/monica/Documents/REU_Materials/EOS_Tables/sly.dat','ap4':'/home/monica/Documents/REU_Materials/EOS_Tables/ap4.dat'}

needs_conversions=['alf1','alf2','sly','ap4']

#read dictionary
def find_eos(model_name):
   if model_name in eos:
      return eos[model_name]
   else:
      print("Please enter valid EOS name")

#import p and rho from table
def p_rho_arrays(model_name):
    file_name = find_eos(model_name)
    dat_file = np.array(np.loadtxt(file_name))

    if model_name in needs_conversions:
      p=dat_file[:,1]*sol**2
      rho=dat_file[:,2] 
    else: 
      p=dat_file[:,2]
      rho=dat_file[:,1]
    
    return p,rho


#iterpolates Log10 of data
def interp_eos_p_of_rho(model_name):
    
    p,rho=p_rho_arrays(model_name)
    n=len(p)
    p=np.log10(p)
    rho=np.log10(rho)

    consts=ms.interpolate(rho,p)
    line_const=ms.lin_extrapolate(rho,p)
    
    #linearly interpolate anything outside range
    line_lower=line_const[0,:]
  
    line_upper=line_const[1,:]

    return consts,line_upper,line_lower

#interpolates Log10 of data
def interp_eos_rho_of_p(model_name):
    
    p,rho=p_rho_arrays(model_name)
    n=len(p)

    p=np.log10(p)
    rho=np.log10(rho)

    consts=ms.interpolate(p,rho)  
    line_const=ms.lin_extrapolate(p,rho)

    #linearly interpolate anything outside range
    line_lower=line_const[0,:]

    line_upper=line_const[1,:]

    return consts,line_upper,line_lower


if __name__ == "__main__":
   print "P(rho) check:"
   
   p,rho=p_rho_arrays('alf1')
   p=np.log10(p)
   rho=np.log10(rho)

   pconst,plineupper,plinelower=interp_eos_p_of_rho('alf1')
  
   for i in np.arange(0,8,1):
      #if pconst[i,0]<2*10**(-16):
      #   pconst[i,0]=0.
      #if pconst[i,1]<2*10**(-16):
      #   pconst[i,1]=0.
      #if pconst[i,2]<2*10**(-16):
      #   pconst[i,2]=0.

      print "cubic"+str(i+1)+"= "+str(p[i])+"+"+str(pconst[i,0])+"(x-"+str(rho[i])+") +"+str(pconst[i,1])+"(x-"+str(rho[i])+")^2 +"+str(pconst[i,2])+"(x-"+str(rho[i])+")^3;"


   print "lineupper="+str(plineupper[0])+"*x+"+str(plineupper[1])
   print "linelower="+str(plinelower[0])+"*x+"+str(plinelower[1])

   print "Piecewise[{",
   for i in np.arange(0,8,1):
       print "{cubic"+str(i+1)+", "+str(rho[i])+"<= x <="+str(rho[i+1])+"}, ",

   print "{linelower, 0<= x <= "+str(rho[0])+"}, {lineupper, "+str(rho[-1])+" <= x <= "+str(rho[-1]+5.)+"}}]"


  # print "rho(P) check:"

   p,rho=p_rho_arrays('alf1')

   p=np.log10(p)
   rho=np.log10(rho)

   rhoconst,rholineupper,rholinelower=interp_eos_rho_of_p('alf1')

   for i in np.arange(0,8,1):
      #if pconst[i,0]<2*10**(-16):
      #   pconst[i,0]=0.
      #if pconst[i,1]<2*10**(-16):
      #   pconst[i,1]=0.
      #if pconst[i,2]<2*10**(-16):
      #   pconst[i,2]=0.

      print "cubic"+str(i+1)+"= "+str(rho[i])+"+"+str(rhoconst[i,0])+"(x-"+str(p[i])+") +"+str(rhoconst[i,1])+"(x-"+str(p[i])+")^2 +"+str(rhoconst[i,2])+"(x-"+str(p[i])+")^3;"

   print "lineupper="+str(rholineupper[0])+"*x+"+str(rholineupper[1])
   print "linelower="+str(rholinelower[0])+"*x+"+str(rholinelower[1])


   print "Piecewise[{",
   for i in np.arange(0,8,1):
       print "{cubic"+str(i+1)+", "+str(p[i])+"<= x <="+str(p[i+1])+"}, ",

   print "{linelower, 0<= x <= "+str(p[0])+"}, {lineupper, "+str(p[-1])+" <= x <= "+str(p[-1]+5.)+"}}]"



