from sys import*
from matplotlib.pyplot import *
matplotlib.use('Agg')
from math import *
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages 
from matplotlib import rc

#rc('text', usetex=True)

# The first part of this files corresponds to functions.
# The main part o the code is indicated by a line:
# ###############################  MAIN CODE  #################################
# Immediatly below, you have a section in which you can define the filter for the selection, etc
# indicated by # PARAMETERS TO SET
# Then, you have a section in which you can change the selection of stars/gal, of the spec-z sample, etc
# indicated by # CONDITION FOR SELECTION



##############################################################
############## DEFINE FUNCTIONS ##############################
##############################################################

#########################
# PLOT THE SPATIAL GALAXY DISTRIBUTION IN SEVERAL REDSHIFT BINS
def area(zmin,zmax):

    #set the axis
    a_min, a_max = np.amin(alpha[cond]), np.amax(alpha[cond])
    d_min, d_max = np.amin(delta[cond]), np.amax(delta[cond])
    axis([a_min,a_max,d_min,d_max])

    # position of all sources
    scatter(alpha[cond], delta[cond], s=1, color='black',alpha=0.5,marker='o')

    # position of the sources within the considered redshift range
    condA = cond & (zs>zmin) & (zs<zmax)
    scatter(alpha[condA], delta[condA], s=3, color='r',marker='^')

    # Labels and titles
    title("$"+str(zmin)+" < z < "+str(zmax)+"$")
    xlabel('$ra$', fontsize=18, labelpad=13)
    ylabel('$dec$', fontsize=18, labelpad=13)

    return


############################
# REDSHIFT DISTRIBUTION OF THE SPEC-Z SAMPLE, ZML AND ZBEST
def distz(magmin,magmax,nstep,ax,leg):

    #set the axis
    ax.axis([0,z_max,0,20.9])

    #new condition with the magnitude range 
    condA = cond  & (mag>magmin) & (mag<magmax) 
    #Check that some object are present
    if len( zs[condA]) > 0 :

       # Histogram with the photometric redshifts minimum chi2
       ax.hist(zs[condA], bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='spectro-z')
     
    # print the legend
    if leg==1:
        ax.legend(prop={'size':8})
    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.01,13),color="black", fontsize=11)

    return


############################
# NUMBER OF FILTERS
def filters(zmin,zmax,nstep,ax) :
    
   #set the axis
   ax.axis([0,19,0,1.9])
   
   #new condition with the magnitude range 
   condA = cond & (zs>zmin) & (zs<zmax)

   nbFil=nbFilt[condA]
   # Histogram with the number of filters 
   ax.hist(nbFil, bins=nstep, histtype='stepfilled', density=1, color='r', label='number of filters')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.1,1),color="black", fontsize=15)
 
   return



############################
# CHI2 DISTRIBUTION
def chi2dist(magmin,magmax,nstep,ax):

    #set the axis
    axis([-2,2.9,0,1.49])

    #new condition with the magnitude range for galaxies, chi gal< chi star
    condA = cond & (mag>magmin) & (mag<magmax) & (nbFilt>2) 
    if len( zs[condA])> 1 :

     chireduit=chi[condA]/(nbFilt[condA]-1)
     # Histogram of log (reduced chi2 gal)
     ax.hist(np.log10(chireduit), bins=nstep, histtype='stepfilled', density=1, color='b')

    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(-1.8,1.2),color="black", fontsize=15)

    return


############################
# MODEL DISTRIBUTION PER REDSHIFT BIN
def model(zmin,zmax,nstep,ax):

    
    #set the axis
    ax.axis([0,19,0,0.99])

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (mod>=0)

    #check if some objects exist
    if len( mod[condA]) > 1 :
        
     # Histogram of the models
     ax.hist(mod[condA], bins=nstep, histtype='stepfilled', density=1, color='b')

    # labels
    ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(1,0.85),color="black", fontsize=15)


    return


############################
# MASS DISTRIBUTION PER REDSHIFT BIN
def massdist(zmin,zmax,nstep,ax,leg):

    
    #set the axis
    ax.axis([6,11.9,0,1.3])

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (massm>0)

    #check if some objects exist
    if len( massm[condA]) > 1 :
        
     # Histogram of the models
     ax.hist(massm[condA], bins=nstep, histtype='stepfilled', density=1, color='b', label='median of the PDF')


    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (massb>0)

    #check if some objects exist
    if len( massb[condA]) > 1 :
        
     # Histogram of the models
     ax.hist(massb[condA], bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='$minimum\; \chi^2$')

    # labels
    ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(7.1,0.85),color="black", fontsize=15)

    # print the legend
    if leg==1:
        ax.legend(prop={'size':8})


    return


############################
# SFR DISTRIBUTION PER REDSHIFT BIN
def sfrdist(zmin,zmax,nstep,ax,leg):

    
    #set the axis
    ax.axis([-4,3,0,0.99])

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (sfrm>-4) & (sfrm<3)

    #check if some objects exist
    if len( sfrm[condA]) > 1 :
        
     # Histogram of the models
     ax.hist(sfrm[condA], bins=nstep, histtype='stepfilled', density=1, color='b', label='median of the PDF')

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (sfrb>-4) & (sfrb<3)

    #check if some objects exist
    if len( sfrb[condA]) > 1 :
        
     # Histogram of the models
     ax.hist(sfrb[condA], bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='$minimum\; \chi^2$')

    # labels
    ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(-3.8,0.85),color="black", fontsize=15)

    # print the legend
    if leg==1:
        ax.legend(prop={'size':8})


    return



############################
# EBV DISTRIBUTION PER REDSHIFT BIN
def ebvDist(zmin,zmax,nstep,ax):

    #set the axis
    ax.axis([0,0.8,0,20])

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (ebv>=0)

    #check if some objects exist
    if len( ebv[condA]) > 1 :
        # Histogram with the E(B-V)
        ax.hist(ebv[condA], bins=nstep, histtype='stepfilled', density=1, color='b')
        # labels
        ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.02,15),color="black", fontsize=15)

    return



############################
# PLOT MASS MEDIAN VERSUS ZBEST
def mass_med_best(magmin,magmax,ax) :
    
   #set the axis
   ax.axis([6,12,6,12])
   
   #new condition with the magnitude range 
   condA = cond & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(massb[condA], massm[condA], s=0.5, color='b',alpha=0.2,marker='s')
   
   #Trace the limits 
   x = np.array([7,12])
   ax.plot(x, x, 'r-')

   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(7.1,10),color="black", fontsize=15)
 
   return


############################
# PLOT SFR MEDIAN VERSUS ZBEST
def sfr_med_best(magmin,magmax,ax) :
    
   #set the axis
   ax.axis([-4,3,-4,3])
   
   #new condition with the magnitude range 
   condA = cond & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(sfrb[condA], sfrm[condA], s=0.5, color='b',alpha=0.2,marker='s')
   
   #Trace the limits 
   x = np.array([-10,12])
   ax.plot(x, x, 'r-')

   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(-3.9,2.5),color="black", fontsize=15)
 
   return


############################
# ABSOLUTE MAGNITUDE IN B BAND VERSUS REDSHIFT
def magabs_z(magmin,magmax,ax):

   #set the axis
   ax.axis([z_min,z_max,-16,-25.9])
   
   #new condition with the magnitude range 
   condA = cond & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(zs[condA], eval("absmag"+str(bFilt))[condA], s=0.5, color='b',alpha=0.2,marker='s')
 
   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.01,-25),color="black", fontsize=15)
 
   return


############################
# MASS VERSUS REDSHIFT
def mass_z(magmin,magmax,ax):

   #set the axis
   ax.axis([z_min,z_max,7,12])
   
   #new condition with the magnitude range 
   condA = cond & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(zs[condA], massm[condA], s=0.5, color='b',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.01,10),color="black", fontsize=15)
 
   return

############################
# SFR VERSUS REDSHIFT
def sfr_z(magmin,magmax,ax):

   #set the axis
   ax.axis([z_min,z_max,-3,4])
   
   #new condition with the magnitude range 
   condA = cond & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(zs[condA], sfrm[condA], s=0.5, color='b',alpha=0.2,marker='s')
 
   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.01,3),color="black", fontsize=15)
 
   return



############################
# SFR VERSUS MASS
def mass_sfr(zmin,zmax,ax):

   #set the axis
   ax.axis([7,12,-3, 4])
   
   #new condition with the magnitude range and star-forming
   condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm>-11) & (ssfrm>-90)
   # Plot photo-z versus spec-z
   ax.scatter(massm[condA],sfrm[condA], s=0.5, color='b',alpha=0.2,marker='s')

   #new condition with the magnitude range and quiescent
   condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm<-11)
   # Plot photo-z versus spec-z
   ax.scatter(massm[condA],sfrm[condA], s=0.5, color='r',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(9,3),color="black", fontsize=15)
 
   return



############################
# lnuv VERSUS SFR
def lnuv_sfr(zmin,zmax,ax):

   #set the axis
   ax.axis([-4, 3,6,12])
   
   #new condition with the magnitude range EBV=0
   condA = cond & (zs>zmin) & (zs<zmax) 
   # be sure to have galaxies in a sufficient number
   if len( massm[condA]) > 3 :
     # Plot sfr versus abs mag
     ax.scatter(sfrm[condA],Lnuv[condA], s=0.5, color='b',alpha=0.2,marker='s')
     
   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(-3,8),color="black", fontsize=15)
 
   return

############################
# ABSOLUTE MAG VERSUS MASS
def absmag_mass(zmin,zmax,ax):

   #set the axis
   ax.axis([7,12,-20,-27])
   
   #new condition with the magnitude range 
   condA = cond & (zs>zmin) & (zs<zmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(massm[condA],eval("absmag"+str(KsFilt))[condA], s=0.5, color='b',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(7.1,-24),color="black", fontsize=15)
 
   return

############################
# ABSOLUTE MAG VERSUS SFR
def absmag_sfr(zmin,zmax,ax):

   #set the axis
   ax.axis([-3,4,-15.1,-25])
   
   #new condition with the magnitude range EBV<0.1
   condA = cond & (zs>zmin) & (zs<zmax) & (ebv<=0.1)   
   # Plot photo-z versus spec-z
   if len( sfrm[condA]) > 0 :
       ax.scatter(sfrm[condA],eval("absmag"+str(uFilt))[condA], s=0.5, color='b',alpha=0.2,marker='s')

   #new condition with the magnitude range 0.1<EBV<0.3
   condA = cond & (zs>zmin) & (zs<zmax) & (ebv>0.1) & (ebv<0.3)   
   # Plot photo-z versus spec-z
   if len( sfrm[condA]) > 0 :
       ax.scatter(sfrm[condA],eval("absmag"+str(uFilt))[condA], s=0.5, color='g',alpha=0.2,marker='s')

   #new condition with the magnitude range 0.3<EBV
   condA = cond & (zs>zmin) & (zs<zmax) & (ebv>=0.3)   
   # Plot photo-z versus spec-z
   if len( sfrm[condA]) > 0 :
       ax.scatter(sfrm[condA],eval("absmag"+str(uFilt))[condA], s=0.5, color='r',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(-3,-21),color="black", fontsize=15)
 
   return


############################
# MASS/LK ratio versus z
def ML_z(magmin,magmax,ax):

   #set the axis
   ax.axis([z_min,z_max,-1.4,0.2])

   # Define the mass-to-light ratio
   MLratio =  massm - Lk + (0.4*(51.605-5.14)) +log10(3.e18*2000./pow(22000.,2.)) - log10(3.826e33)

   #new condition with the magnitude range and star-forming
   condA = cond & (mag>magmin) & (mag<magmax) & (ssfrm>-11) & (ssfrm>-90)
   if len( zs[condA]) > 0 :
     # Plot photo-z versus spec-z
     ax.scatter(zs[condA],MLratio[condA], s=0.5, color='b',alpha=0.2,marker='s')

   #new condition with the magnitude range and quiescent
   condA = cond & (mag>magmin) & (mag<magmax) & (ssfrm<-11)
   if len( zs[condA]) > 0 :
     # Plot photo-z versus spec-z
     ax.scatter(zs[condA],MLratio[condA], s=0.5, color='r',alpha=0.2,marker='s')


   #Trace the limits 
   x = np.array([0,6])
   ax.plot(x, -0.27*x-0.05 - 0.24, 'r-')
   ax.plot(x, -0.18*x-0.05 - 0.24, 'r-')
     
   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.01,0.1),color="black", fontsize=15)
 
   return


############################
# COLOR VERUS M PLOT 
def RFcolor(zmin,zmax,ax) :
    
   #set the axis
   ax.axis([-25,-12,-1,3])
   

   # check that the U band and R band filters are well defined
   if ((uFilt>=0) & (rFilt>=0)):
    magx= eval("absmag"+str(rFilt))
    coly= eval("absmag"+str(uFilt)+"-"+"absmag"+str(rFilt))
    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm>-11) & (ssfrm>-90)
    # Scatter plot 
    ax.scatter(magx[condA], coly[condA], s=0.5, color='b',alpha=0.2,marker='s')

    #new condition with the magnitude range 
    condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm<-11) & (ssfrm>-90)
    # Scatter plot 
    ax.scatter(magx[condA], coly[condA], s=0.5, color='r',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(-24,2),color="black", fontsize=15)
 
   return

############################
# COLOR-COLOR PLOT A LA WILLIAMS
def william(zmin,zmax,ax) :
    
   #set the axis
   ax.axis([-0.1,2.4,-1,3])
   

   # check that the U, R, J bands exists
   if ((uFilt>=0) &  (rFilt>=0) & (jFilt>=0)) :
    colx= eval("absmag"+str(rFilt)+"-"+"absmag"+str(jFilt))
    coly= eval("absmag"+str(uFilt)+"-"+"absmag"+str(rFilt))
    
    #new condition with the magnitude range and star-forming
    condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm>-11) & (ssfrm>-90)
    # Scatter plot color-color
    ax.scatter(colx[condA], coly[condA], s=0.5, color='b',alpha=0.2,marker='s')

    #new condition with the magnitude range and quiescent
    condA = cond & (zs>zmin) & (zs<zmax) & (ssfrm<-11)  & (ssfrm>-90)
    # Scatter plot color-color
    ax.scatter(colx[condA], coly[condA], s=0.5, color='r',alpha=0.2,marker='s')
    
   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.1,2),color="black", fontsize=15)
 
   return


############################
# CHECK ERROR AT 68%
def check_error() :
    
   #set the axis
   axis([15,27,-0.4,0.4])

   # 1 sigma error on the photo-z
   diffu=masss-massm
   diffl=massl-massm
   mag=eval("mag"+str(selFilt))

   condA = cond
   
   # Error versus mag
   scatter(mag[condA], diffu[condA], s=0.5, color='b',alpha=0.2,marker='s')
   scatter(mag[condA], diffl[condA], s=0.5, color='r',alpha=0.2,marker='s')

   # 68%
   plot([-10,100], [0,0.], 'r--')

   xlabel('magnitude', fontsize=18, labelpad=13)
   ylabel('$68\% \; mass \; error$', fontsize=18, labelpad=13)

   return


############################
# ERROR ON THE PHOTOZ VERSUS MAG
def errorMag(zmin,zmax,ax):

   # step en mag
   mstep=0.2
   mstart=16
   mnbstep=150

   medpVec=[]
   mednVec=[]
   magVec=[]
   
   #set the axis
   ax.axis([mag_min,mag_max,-0.4,0.4])

   # Loop over the magnitude bins
   for rm in range(mnbstep) :

      # take the maximum of two vector
      diffp=masss-massm
      diffn=massl-massm

      #limit inf and sup
      mlimi=mstart+mstep*rm
      mlims=mstart+mstep*(rm+1)

      #new condition with the magnitude range 
      condA = cond & (mag>mlimi) & (mag<=mlims) & (zs>zmin) & (zs<zmax) & (massm>0)

      # be sure to have galaxies in a sufficient number
      if len( zs[condA]) > 3 :
       
       #compute the median value
       medp=np.median(diffp[condA])
       medn=np.median(diffn[condA])
       
       # Add the result to the vector
       medpVec.append(medp)
       mednVec.append(medn)
       magVec.append(mstart+mstep*(rm+0.5))
      
   # Histogram with the photometric redshifts zml
   ax.plot(magVec,medpVec, color='b',linestyle='-')
   ax.plot(magVec,mednVec, color='r',linestyle='-')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(mstart,0.60),color="black", fontsize=15)


############################
# ERROR ON THE PHOTOZ VERSUS REDSHIFT
def errorZ(magmin,magmax,ax):

   # step en mag
   zstep=0.02
   zstart=0
   znbstep=60

   medpVec=[]
   mednVec=[]
   zVec=[]
   
   #set the axis
   ax.axis([z_min,z_max,-0.4,0.4])

   # Loop over the magnitude bins
   for rm in range(znbstep) :

      # take the maximum of two vector
      diffp=masss-massm
      diffn=massl-massm

      #limit inf and sup
      zlimi=zstart+zstep*rm
      zlims=zstart+zstep*(rm+1)

      #new condition with the magnitude range 
      condA = cond & (zs>zlimi) & (zs<=zlims) & (mag>magmin) & (mag<magmax) & (massm>0)

      # be sure to have galaxies in a sufficient number
      if len( zs[condA]) > 3 :
       
       #compute the median value
       medp=np.median(diffp[condA])
       medn=np.median(diffn[condA])
       
       # Add the result to the vector
       medpVec.append(medp)
       mednVec.append(medn)
       zVec.append(zstart+zstep*(rm+0.5))
      
   # Histogram with the photometric redshifts zml
   ax.plot(zVec,medpVec, color='b',linestyle='-')
   ax.plot(zVec,mednVec, color='r',linestyle='-')

   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(zstart,0.60),color="black", fontsize=15)

   
############################
# ERROR ON THE MAG VERSUS MAG
def errorMagMag(flt,ax,labFl):

   
   #set the axis
   ax.axis([20,29,0,0.5])

   mag=eval("mag"+str(flt)+"[cond]")
   emag=eval("emag"+str(flt)+"[cond]")
   ax.scatter(mag,emag,s=0.5, color='b',alpha=0.2,marker='s')


   # labels
   ax.annotate(labFl,xy=(21,0.20),color="black", fontsize=15)


       
#############################################################################
###############################  MAIN CODE  #################################
#############################################################################

########
# PARAMETERS TO SET

# Number of the filter start at 0
selFilt=5   # filter for the selection in mag
uFilt=0
bFilt=1
rFilt=3
zFilt=5
jFilt=8
KsFilt=10

# Array in redshift and mag, isolate extreme values
range_z   = [0,1,2,3,6]
z_min, z_max = np.amin(range_z), np.amax(range_z)
range_mag = [15.,22.5,23.5,25,28]
mag_min, mag_max = np.amin(range_mag), np.amax(range_mag)


# define the sub panels
# case with several bin of mag
if(len(range_mag)>1):
    # Number of row with 2 columns
    nbRowM=int(ceil(float(len(range_mag)-1)/2.))
    print("Number of lines ",len(range_mag),(float(len(range_mag))/2.),(ceil(float(len(range_mag))/2.)),nbRowM)
    nbColM=2
else:
    nbColM=1
    nbRowM=1

# define the sub panels
# case with several bin of z
if(len(range_z)>1):
    # Number of row with 2 columns
    nbRowZ=int(ceil(float(len(range_z)-1)/2.))
    nbColZ=2
else:
    nbColZ=1
    nbRowZ=1
    
# Fraction of sources to select
fracSub=0.1



########
# READ THE INPUT FILE

# Read the first argument with the name of the Le Phare output catalogue
fileIn=sys.argv[1]
catIn=open(fileIn,'r')
print("Name of the Le Phare output catalogue : ",fileIn)


# Loop over the filters
nbFilt=30
magst=""
idmagst=""
# create the string to read the absolute mag
for i in range(nbFilt) :
    magst=magst+",absmag"+str(i)
    idmagst=idmagst+","+str(i+30)
# create the string to read the mag
for i in range(nbFilt) :
    magst=magst+",mag"+str(i)
    idmagst=idmagst+","+str(i+30+nbFilt)
# create the string to read the error mag
for i in range(nbFilt) :
    magst=magst+",emag"+str(i)
    idmagst=idmagst+","+str(i+30+2*nbFilt)


# Extract from the ascii file
commandst = "Id,zs,mod,ebv,law,chi,nbFilt,ageb,agel,agem,ages,ldustb,ldustl,ldustm,ldusts,massb,massl,massm,masss,sfrb,sfrl,sfrm,sfrs,ssfrb,ssfrl,ssfrm,ssfrs,Lnuv,Lr,Lk"+magst+"= np.loadtxt(catIn, dtype='float', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29"+idmagst+"), unpack=True )"



print(commandst)

# transform the string into a command
exec(commandst)


##################################################
#    CONDITION FOR SELECTION

# Mag use to select the sample
mag=eval("mag"+str(selFilt))

# Generate a random number to subselect a sample
subsample = np.random.uniform(low=0., high=1., size=(len( zs),))

# General condition to select the galaxies in the expected z/mag range
#cond = (zs>z_min) & (zs<z_max) & (mag>mag_min) & (mag<mag_max) & (mask==0) & (subsample<fracSub) & (fl>-0.1) & (fl<0.1) & (nbFilt>-1)
cond = (zs>z_min) & (zs<z_max) & (mag>mag_min) & (mag<mag_max)  & (nbFilt>-1)

 
##################################################
#### START THE FIGURES

# All the figures will be collected in a single pdf file 
pdfOut = PdfPages('figuresLPP.pdf')


##################
######## AREA

           
# Loop over the magnitude bins
#for rm in range(len(range_z)) :

   #if (rm>0):
    # Define the subplots and pass the panel
    #clf()
    #area(range_z[rm-1],range_z[rm])

    # store the figure in a PDF
    #savefig(pdfOut,format='pdf')
   

#############################
####### REDSHIFT DISTRIBUTION


# redshift step for the histogram
nstep=20

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f,axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$Redshift$', ha='center')
f.text(0.04, 0.5, '$N_{normalized}$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   # No legend
   leg=0
   
   if (rm>0):
    if rm==1:
        leg=1  #legend in the first panel
    # Define the subplots and pass the panel
    distz(range_mag[rm-1],range_mag[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)],leg)

# store the figure in a PDF
savefig(pdfOut,format='pdf')




#############################
######## NUMBER OF FILTERS USED IN THE FIT
#
#print 'number of filters'
#
## redshift step for the histogram
#nstep=10
#
### Create a figure with an array with nbRowM*nbColM subpanels
#clf()
#f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
## No space between the figures
#subplots_adjust( hspace=0,wspace=0)
#
## label of the figure
#f.text(0.5, 0.04, 'Number of filters', ha='center')
#f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
#                     
## Loop over the redshift bins
#for rm in range(len(range_z)) :
#
#   if (rm>0):
#    # Define the subplots and pass the panel
#    filters(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])
#
## store the figure in a PDF
#savefig(pdfOut,format='pdf')


############################
##### CHI2 DISTRIBUTION PER MAGNITUDE BIN



# nb ebv step for the histogram
nstep=20

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log($\chi^2$)', ha='center')
f.text(0.04, 0.5, '$N_{normalized}$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   # No legend
   leg=0
   
   if (rm>0):
    if rm==2:
        leg=1  #legend in the first panel
    # Define the subplots and pass the panel
    chi2dist(range_mag[rm-1],range_mag[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### MODEL


# redshift step for the histogram
nstep=20

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'best fit template', ha='center')
f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    model(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### MASS


# mass step for the histogram
nstep=10

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(stellar mass)', ha='center')
f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   # No legend
   leg=0

   if (rm>0):
    if rm==1:
        leg=1  #legend in the first panel

    # Define the subplots and pass the panel
    massdist(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)],leg)

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### SFR


# mass step for the histogram
nstep=10

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(SFR)', ha='center')
f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   # No legend
   leg=0

   if (rm>0):
    if rm==4:
        leg=1  #legend in the first panel

    # Define the subplots and pass the panel
    sfrdist(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)],leg)

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### EBV DISTRIBUTION


# number of steps for the histogram
nstep=10

## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'E(B-V)', ha='center')
f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    ebvDist(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')



############################
# PLOT MASS MEDIAN VERSUS ZBEST



## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$log(mass) \; minimum\; \chi^2$', ha='center')
f.text(0.04, 0.5, '$log(mass)\; median\; PDZ$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    mass_med_best(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])


# store the figure in a PDF
savefig(pdfOut,format='pdf')

############################
# PLOT SFR MEDIAN VERSUS ZBEST



## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$log(SFR) \; minimum\; \chi^2$', ha='center')
f.text(0.04, 0.5, '$log(SFR)\; median\; PDZ$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    sfr_med_best(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])


# store the figure in a PDF
savefig(pdfOut,format='pdf')
   
#############################
####### sfr versus mass


## Create a figure with an array with nbRowZ*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(mass)', ha='center')
f.text(0.04, 0.5, 'log(sfr)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    mass_sfr(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### mag abs versus mass


## Create a figure with an array with nbRowZ*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(mass)', ha='center')
f.text(0.04, 0.5, 'log(absolute magnitude K)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    absmag_mass(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### ML versus z


## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(z)', ha='center')
f.text(0.04, 0.5, 'log(mass-to-light ratio)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    ML_z(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### mag abs versus sfr


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(SFR)', ha='center')
f.text(0.04, 0.5, 'log(absolute magnitude U)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    absmag_sfr(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')



#############################
####### mag abs versus sfr


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'log(SFR)', ha='center')
f.text(0.04, 0.5, 'log(L NUV corrected for extinction in erg/s/Hz)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    lnuv_sfr(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')



#############################
####### absolute magnitude versus redshift


## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$Redshift$', ha='center')
f.text(0.04, 0.5, '$M_B$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    magabs_z(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')

#############################
####### mass versus redshift


## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$Redshift$', ha='center')
f.text(0.04, 0.5, 'log(mass)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    mass_z(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')

#############################
####### sfr versus redshift


## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$Redshift$', ha='center')
f.text(0.04, 0.5, 'log(sfr)', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    sfr_z(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


###############################
######### REST-FRAME COLORS VERSUS MAG



## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$M_R$', ha='center')
f.text(0.04, 0.5, '$M_{U}-M_R$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    RFcolor(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')

##############################
######## WILLIAM PLOT


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$M_R-M_J$', ha='center')
f.text(0.04, 0.5, '$M_{U}-M_R$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    william(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')



############################
##### ERROR 68% versus mag


# clear the figure
clf()
check_error()
# store the figure in a PDF
savefig(pdfOut,format='pdf')


###############################
######### MASS ERRORS VERSUS I


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$mag$', ha='center')
f.text(0.04, 0.5, '$mass \; uncertainties$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    errorMag(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


###############################
######### MASS ERRORS VERSUS Z


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$z$', ha='center')
f.text(0.04, 0.5, '$mass \; uncertainties$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    errorZ(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')

################################
########## MAG ERRORS VERSUS MAG
#
#
### Create a figure with an array with nbRowZ*nbColZ subpanels
#clf()
#f, axarr = subplots(2,2,sharex=True,sharey=True)
## No space between the figures
#subplots_adjust( hspace=0,wspace=0)
#
## label of the figure
#f.text(0.5, 0.04, '$mag$', ha='center')
#f.text(0.04, 0.5, '$mag \; uncertainties$', va='center', rotation='vertical')
#
#
## Define the subplots and pass the panel
#errorMagMag(0,axarr[0,0],'125W')
#errorMagMag(1,axarr[1,0],'140W')
#errorMagMag(2,axarr[0,1],'160W')
#errorMagMag(6,axarr[1,1],'606W')
#
## store the figure in a PDF
#savefig(pdfOut,format='pdf')

pdfOut.close()


