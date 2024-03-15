from sys import*
from matplotlib.pyplot import *
matplotlib.use('Agg')
from math import *
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages 
from matplotlib import rc


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
    condA = cond & (zp>zmin) & (zp<zmax)
    scatter(alpha[condA], delta[condA], s=3, color='r',marker='^')

    # Labels and titles
    title("$"+str(zmin)+" < z < "+str(zmax)+"$")
    xlabel('$ra$', fontsize=18, labelpad=13)
    ylabel('$dec$', fontsize=18, labelpad=13)

    return


############################
# PLOT PHOTO-Z VERSUS SPEC-Z
def zp_zs(magmin,magmax) :
    
   #set the axis
   axis([z_min,z_max,z_min,z_max])

   #new condition with the magnitude range 
   condA = cond & condspec & (mag>magmin) & (mag<magmax)  

   # Plot photo-z versus spec-z
   scatter(zs[condA], zp[condA], s=1, color='b',alpha=0.5,marker='s')

   #Check that we have some sources before performing statistics
   Ngal = len( zp[condA])
   if Ngal > 0 :

       #statistics
       arg_bias = zp[condA]-zs[condA]
       arg_std = arg_bias / (1. + zs[condA])
       NMAD = 1.4821 * np.median( abs(arg_std))
       cond_outl = ( abs(arg_std) > 0.15 )
       outl_rate = len(arg_std[cond_outl]) / float(Ngal)
       annotate(r'$N_{gal}  = '+str(Ngal)+'$ \n'+'$\eta  ='+str(100*round(outl_rate,3))+'  \%$\n'+'$ \sigma_{\Delta z /(1+z)}  = '+str(round(NMAD,5))+'$',xy=(0.1*z_max,0.8*z_max),color="black", fontsize=15)

   #Trace the limits 0.15(1+z)
   x_zs = np.array([0,6])
   plot(x_zs, x_zs*1.15+0.15, 'c--')
   plot(x_zs, x_zs, 'r-')
   plot(x_zs, x_zs*0.85-0.15, 'c--')

   #labels
   title("$"+str(magmin)+" < mag < "+str(magmax)+"$")
   xlabel('$z_{spec}$', fontsize=18, labelpad=13)
   ylabel('$z_{phot} \; minimum \; \chi^2$', fontsize=18, labelpad=13)

   return


############################
# PLOT PHOTO-Z ZML VERSUS ZBEST
def zml_zp(magmin,magmax,ax) :
    
   #set the axis
   ax.axis([z_min,z_max,z_min,z_max])
   
   #new condition with the magnitude range 
   condA = cond & condgal & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(zp[condA], zml[condA], s=1, color='b',alpha=0.5,marker='s')
   
   #Trace the limits 0.15(1+z)
   x_zs = np.array([0,6])
   ax.plot(x_zs, x_zs*1.15+0.15, 'c--')
   ax.plot(x_zs, x_zs, 'r-')
   ax.plot(x_zs, x_zs*0.85-0.15, 'c--')

   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.1*z_max,0.8*z_max),color="black", fontsize=15)
 
   return



############################
# REDSHIFT DISTRIBUTION OF THE SPEC-Z SAMPLE, ZML AND ZBEST
def distz(magmin,magmax,nstep,ax,leg):

    #set the axis
    ax.axis([0,z_max,0,5.9])

    #new condition with the magnitude range 
    condA = cond  & condgal & (mag>magmin) & (mag<magmax) & (zml>0)
    #Check that some object are present
    if len( zp[condA]) > 0 :

       # Histogram with the photometric redshifts median PDF
       ax.hist(zml[condA], bins=nstep, histtype='stepfilled', density=1, color='b', label='$z_{phot}\; median\; PDZ$')

       # Histogram with the photometric redshifts minimum chi2
       ax.hist(zp[condA], bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='$z_{phot}\; minimum\; \chi^2$')

     
    #new condition with the magnitude range and spectro-z
    condA = cond & condspec & (mag>magmin) & (mag<magmax) 
    #Check that some object are present
    if len( zp[condA]) > 0 :

        # Histogram with the photometric redshifts zbest
        ax.hist(zs[condA], bins=nstep, histtype='stepfilled', density=1, color='g', alpha=0.2, label='$z_{spectro}$')

    # print the legend
    if leg==1:
        ax.legend(prop={'size':8})
    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.05*z_max,4.5),color="black", fontsize=11)

    return


############################
# REDSHIFT DISTRIBUTION OF THE SECOND PEAK PHOTO-Z, ZML AND ZBEST OF THESE SOURCES
def secondpeak(magmin,magmax,nstep,ax,leg):

    #set the axis
    axis([0,z_max,0,5])

    #new condition with the magnitude range 
    condA = cond & condgal & (mag>magmin) & (mag<magmax) & (zp2>0)

    # If some sources
    if len( zp2[condA]) > 1 :

     # Histogram with the photometric redshifts median PDF
     ax.hist(zp2[condA], bins=nstep, histtype='stepfilled', density=1, color='b', label='$z_{phot}\; second\; peak$')

     # Histogram with the photometric redshifts minimum chi2
     ax.hist(zp[condA], bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='$z_{phot}\; minimum \; \chi^2$')

    # print the legend
    if leg==1:
        ax.legend(prop={'size':8})

    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.2,3.5),color="black", fontsize=13)

    return


############################
# NUMBER OF FILTERS
def filters(zmin,zmax,nstep,ax) :
    
   #set the axis
   ax.axis([0,50,0,3])
   
   #new condition with the magnitude range 
   condA = cond & (zp>zmin) & (zp<zmax)

   # Histogram with the number of filters 
   ax.hist(nbFilt, bins=nstep, histtype='stepfilled', density=1, color='r')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.1,1),color="black", fontsize=15)
 
   return



############################
# CHI2 DISTRIBUTION
def chi2dist(magmin,magmax,nstep,ax,leg):

    #set the axis
    axis([-2,2.9,0,1.49])

    #new condition with the magnitude range for galaxies, chi gal< chi star
    condA = cond & condgal & (mag>magmin) & (mag<magmax) & (nbFilt>2) & (chi>0)
    if len( zp[condA])> 1 :

     chireduit=chi[condA]/(nbFilt[condA]-1)
     # Histogram of log (reduced chi2 gal)
     ax.hist(np.log10(chireduit), bins=nstep, histtype='stepfilled', density=1, color='b', label='$\chi^2_{gal} \; if \; \chi^2_{gal}<\chi^2_{stars}$')

    #new condition with the magnitude range for stars, chi star < chi gal
    condA = cond & condstar & (mag>magmin) & (mag<magmax) & (nbFilt>2)  & (chi>0)
    if len( zp[condA]) > 1 :

     chireduit=chis[condA]/(nbFilt[condA]-1)
     # Histogram of log (reduced chi2 gal)
     ax.hist(np.log10(chireduit), bins=nstep, histtype='stepfilled', density=1, color='r', alpha=0.5, label='$\chi^2_{star} \; if \;  \chi^2_{gal}>\chi^2_{stars}$')

    # print the legend
    if leg==1:
      ax.legend(prop={'size':8})
    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(-1.8,1.2),color="black", fontsize=15)

    return


############################
# MODEL DISTRIBUTION PER REDSHIFT BIN
def model(zmin,zmax,nstep,ax):

    
    #set the axis
    ax.axis([0,70,0,0.99])

    #new condition with the magnitude range 
    condA = cond & condgal & (zp>zmin) & (zp<zmax)

    #check if some objects exist
    if len( mod[condA]) > 1 :

     # Histogram of the models
     ax.hist(mod[condA], bins=nstep, histtype='stepfilled', density=1, color='b')

     # labels
     ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(1,0.15),color="black", fontsize=15)


    return



############################
# EBV DISTRIBUTION PER REDSHIFT BIN
def ebvDist(zmin,zmax,nstep,ax):

    #set the axis
    ax.axis([0,0.6,0,20])

    #new condition with the magnitude range 
    condA = cond & condgal & (zp>zmin) & (zp<zmax)

    #check if some objects exist
    if len( ebv[condA]) > 1 :

     # Histogram with the E(B-V)
     ax.hist(ebv[condA], bins=nstep, histtype='stepfilled', density=1, color='b')

     # labels
     ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.02,15),color="black", fontsize=15)

    return


############################
# BzK
def BzK(zmin,zmax,ax,galstar):

    #set the axis
    ax.axis([-0.5, 5.5, -1., 4.5])

    #new condition with the magnitude range
    if galstar==0:
      condA = cond & condgal & (zp>zmin) & (zp<zmax)  
    else:
      condA = cond & condstar

    if ((bFilt>=0) & (zFilt>=0) & (KsFilt>=0)):
      colx= eval("mag"+str(bFilt)+"-"+"mag"+str(zFilt))
      coly= eval("mag"+str(zFilt)+"-"+"mag"+str(KsFilt))
      # Scatter plot BzK
      ax.scatter(colx[condA], coly[condA], s=0.5, color='b',alpha=0.2,marker='s')

    # BzK limit
    ax.plot([-1,6], [-0.5,6])

    # labels
    if galstar==0:
      ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.1,3),color="black", fontsize=15)
    else:
      ax.annotate("$\chi_s<\chi_g$",xy=(0.1,3),color="black", fontsize=15)

    return



############################
# ABSOLUTE MAGNITUDE IN B BAND VERSUS REDSHIFT
def magabs_z(magmin,magmax,ax):

   #set the axis
   ax.axis([z_min,z_max,-16,-25.9])
   
   #new condition with the magnitude range 
   condA = cond & condgal & (mag>magmin) & (mag<magmax)
   
   # Plot photo-z versus spec-z
   ax.scatter(zp[condA], eval("absmag"+str(bFilt))[condA], s=0.5, color='b',alpha=0.2,marker='s')
 
   # labels
   ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.1,-25),color="black", fontsize=15)
 
   return

############################
# COLOR VERUS M PLOT 
def RFcolor(zmin,zmax,ax) :
    
   #set the axis
   ax.axis([-24.9,-16,-0.7,2.4])
   
   #new condition with the magnitude range 
   condA = cond & condgal & (zp>zmin) & (zp<zmax)

   # check that the U band and R band filters are well defined
   if ((uFilt>=0) & (rFilt>=0)):
    magx= eval("absmag"+str(rFilt))
    coly= eval("absmag"+str(uFilt)+"-"+"absmag"+str(rFilt))
    
    # Scatter plot 
    ax.scatter(magx[condA], coly[condA], s=0.5, color='b',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(-20,1),color="black", fontsize=15)
 
   return

############################
# COLOR-COLOR PLOT A LA WILLIAMS
def william(zmin,zmax,ax) :
    
   #set the axis
   ax.axis([-2.1,1.9,0,2.5])
   
   #new condition with the magnitude range 
   condA = cond & condgal & (zp>zmin) & (zp<zmax) 

   # check that the U, R, J bands exists
   if ((uFilt>=0) &  (rFilt>=0) & (jFilt>=0)) :
    colx= eval("absmag"+str(rFilt)+"-"+"absmag"+str(jFilt))
    coly= eval("absmag"+str(uFilt)+"-"+"absmag"+str(rFilt))
    
    # Scatter plot color-color
    ax.scatter(colx[condA], coly[condA], s=0.5, color='b',alpha=0.2,marker='s')

   # labels
   ax.annotate("$"+str(zmin)+" < z < "+str(zmax)+"$",xy=(0.1,2),color="black", fontsize=15)
 
   return


############################
# CUMULATIVE DISTRIBUTION ERROR ZP
def cumulative68(magmin,magmax,nstep,ax):

    #set the axis
    ax.axis([0,9.9,0,0.99])
    
    # take the maximum of two vectors in the mdedian PDF case
    diffzml=abs(zml-zs)/np.maximum(abs(zmlu68-zml),abs(zmll68-zml))
    # take the maximum of two vectors in the minimum chi2 case
    diffzp=abs(zp-zs)/np.maximum(abs(zu68-zp),abs(zl68-zp))

    #new condition with the magnitude range 
    condAzml = cond & condspec & (mag>magmin) & (mag<magmax) & (diffzml<10) & (zml>0)
    condAzp = cond & condspec & (mag>magmin) & (mag<magmax) & (diffzp<10) & (zp>0)

    # some galaxy exists
    if len( zml[condAzml]) > 1 :
    
       # Histogram with the photometric redshifts zml
       ax.hist(diffzp[condAzp], bins=nstep, histtype='step', density=1, color='b', cumulative=True, label='$z_{phot}\; minimum \; \chi^2$')
       # Histogram with the photometric redshifts zp
       ax.hist(diffzml[condAzml], bins=nstep, histtype='step', density=1, color='r', cumulative=True, label='$z_{phot}\; median \; PDF$')

       # 68%
       ax.plot([1,1,0], [0,0.68,0.68], color='r')

       # legend
       ax.legend(prop={'size':8})

    # labels
    ax.annotate("$"+str(magmin)+" < mag < "+str(magmax)+"$",xy=(0.1,0.80),color="black", fontsize=15)

    
    return

############################
# CHECK ERROR AT 68%
def check_error() :
    
   #set the axis
   axis([15,27,-1,1])

   # 1 sigma error on the photo-z
   diffu=zu68-zp
   diffl=zl68-zp
   mag=eval("mag"+str(selFilt))

   condA = cond & condgal
   
   # Error versus mag
   scatter(mag[condA], diffu[condA], s=1, color='b',alpha=0.5,marker='s')
   scatter(mag[condA], diffl[condA], s=1, color='r',alpha=0.5,marker='s')

   # 68%
   plot([-10,100], [0,0.], color='r')

   xlabel('magnitude', fontsize=18, labelpad=13)
   ylabel('$68\% \; photo-z \; error$', fontsize=18, labelpad=13)

   return


############################
# ERROR ON THE PHOTOZ VERSUS MAG
def errorMag(zmin,zmax,ax):

   # step en mag
   mstep=0.2
   mstart=16
   mnbstep=100

   medpVec=[]
   mednVec=[]
   magVec=[]
   
   #set the axis
   ax.axis([15,27,-1,1])

   # Loop over the magnitude bins
   for rm in range(mnbstep) :

      # take the maximum of two vector
      diffp=zmlu68-zml
      diffn=zmll68-zml

      #limit inf and sup
      mlimi=mstart+mstep*rm
      mlims=mstart+mstep*(rm+1)

      #new condition with the magnitude range 
      condA = cond & condgal & (mag>mlimi) & (mag<=mlims) & (zml>zmin) & (zml<zmax)

      # be sure to have galaxies in a sufficient number
      if len( zml[condA]) > 3 :
       
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
   zstep=0.1
   zstart=0
   znbstep=60

   medpVec=[]
   mednVec=[]
   zVec=[]
   
   #set the axis
   ax.axis([0,5.9,-1,1])

   # Loop over the magnitude bins
   for rm in range(znbstep) :

      # take the maximum of two vector
      diffp=zmlu68-zml
      diffn=zmll68-zml

      #limit inf and sup
      zlimi=zstart+zstep*rm
      zlims=zstart+zstep*(rm+1)

      #new condition with the magnitude range 
      condA = cond & condgal & (zml>zlimi) & (zml<=zlims) & (mag>magmin) & (mag<magmax)

      # be sure to have galaxies in a sufficient number
      if len( zml[condA]) > 3 :
       
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

       
#############################################################################
###############################  MAIN CODE  #################################
#############################################################################

########
# PARAMETERS TO SET

# Number of the filter start at 0
selFilt= 4  # filter for the selection in mag
uFilt=0
bFilt=2
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

########
# READ THE INPUT FILE

# Read the first argument with the name of the photo-z output catalogue
fileIn=sys.argv[1]
catIn=open(fileIn,'r')
print("Name of the photo-z catalogue : ",fileIn)


# Loop over the filters
nbFilt=30
magst=""
idmagst=""
# create the string to read the mag
for i in range(nbFilt) :
    magst=magst+",mag"+str(i)
    idmagst=idmagst+","+str(i+20)
# create the string to read the error mag
for i in range(nbFilt) :
    magst=magst+",emag"+str(i)
    idmagst=idmagst+","+str(i+20+nbFilt)
# create the string to read the absolute mag
for i in range(nbFilt) :
    magst=magst+",absmag"+str(i)
    idmagst=idmagst+","+str(i+20+3*nbFilt)
# create the string to read the uncertainties on absolute mag
for i in range(nbFilt) :
    magst=magst+",eabsmag"+str(i)
    idmagst=idmagst+","+str(i+20+4*nbFilt)


# Extract from the ascii file
commandst = "Id,zp,zl68,zu68,zml,zmll68,zmlu68,chi,mod,law,ebv,zp2,chi2,mod2,ebv2,zq,chiq,modq,mods,chis"+magst+",scale,nbFilt,context,zs= np.loadtxt(catIn, dtype='float', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19"+idmagst+",200,201,202,203), unpack=True )"
# transform the string into a command
print(commandst)
exec(commandst)


##################################################
#    CONDITION FOR SELECTION

# Mag use to select the sample
mag=eval("mag"+str(selFilt))

# General condition to select the galaxies in the expected z/mag range
cond = (zp>z_min) & (zp<z_max) & (mag>mag_min) & (mag<mag_max) 

# condition to select stars
condstar = (chis<chi)

# condition to select galaxies
condgal =  (~condstar)

# condition to select spectroscopic redshifts
condspec = (zs>0) & (zs<9)

 
##################################################
#### START THE FIGURES

# All the figures will be collected in a single pdf file 
pdfOut = PdfPages('figuresLPZ.pdf')


   
############################
##### PHOTO-Z VERSUS SPEC-Z 


# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   # define the variable range_min, max
   if (rm>0):
    # clear the figure
    clf()
    zp_zs(range_mag[rm-1],range_mag[rm])

    # store the figure in a PDF
    savefig(pdfOut,format='pdf', bbox_inches='tight')


############################
##### PHOTO-Z ZML VERSUS ZBEST


## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$z_{phot}\; minimum\; \chi^2$', ha='center')
f.text(0.04, 0.5, '$z_{phot}\; median\; PDZ$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    zml_zp(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])


# store the figure in a PDF
savefig(pdfOut,format='pdf')
   

#############################
####### REDSHIFT DISTRIBUTION


# redshift step for the histogram
nstep=20

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
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




#########################################
####### REDSHIFT DISTRIBUTION SECOND PEAK


# redshift step for the histogram
nstep=20

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
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
    if rm==4:
        leg=1  #legend in the first panel

    # Define the subplots and pass the panel
    secondpeak(range_mag[rm-1],range_mag[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)],leg)

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
####### NUMBER OF FILTERS USED IN THE FIT


# redshift step for the histogram
nstep=10

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, 'Number of filters', ha='center')
f.text(0.04, 0.5, 'N', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    filters(range_z[rm-1],range_z[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


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
    chi2dist(range_mag[rm-1],range_mag[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)],leg)

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
##### BzK per redshift bin



## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$B-z$', ha='center')
f.text(0.04, 0.5, '$z-K$', va='center', rotation='vertical')
                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    BzK(0.001,1.4,axarr[0,0],0)
    BzK(1.4,2.7,axarr[0,1],0)
    BzK(2.7,6,axarr[1,0],0)
    BzK(-0.001,6,axarr[1,1],1)


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


###############################
######### REST-FRAME COLORS VERSUS MAG



## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$M_R$', ha='center')
f.text(0.04, 0.5, '$M_U-M_R$', va='center', rotation='vertical')
                     
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
f.text(0.04, 0.5, '$M_U-M_R$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    william(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


#############################
###### CUMULATIVE 68


# nb ebv step for the histogram
nstep=40

## Create a figure with an array with nbRowM*nbColM subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$abs((z_{spec}-z_{phot})/z_{phot} uncertainties)$', ha='center')
f.text(0.04, 0.5, '$cumulative(N)$', va='center', rotation='vertical')

                     
# Loop over the magnitude bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    cumulative68(range_mag[rm-1],range_mag[rm],nstep,axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

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
######### PHOTO-Z ERRORS VERSUS I


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowZ,nbColZ,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$mag$', ha='center')
f.text(0.04, 0.5, '$z_{phot} \; uncertainties$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_z)) :

   if (rm>0):
    # Define the subplots and pass the panel
    errorMag(range_z[rm-1],range_z[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')


###############################
######### PHOTO-Z ERRORS VERSUS Z


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(nbRowM,nbColM,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$z_{phot}$', ha='center')
f.text(0.04, 0.5, '$z_{phot} \; uncertainties$', va='center', rotation='vertical')
                     
# Loop over the redshift bins
for rm in range(len(range_mag)) :

   if (rm>0):
    # Define the subplots and pass the panel
    errorZ(range_mag[rm-1],range_mag[rm],axarr[int(ceil(rm/2.)-1),int(ceil(rm % 2)-1)])

# store the figure in a PDF
savefig(pdfOut,format='pdf')

pdfOut.close()


