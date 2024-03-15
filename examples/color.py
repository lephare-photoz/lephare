from sys import*
from matplotlib.pyplot import *
from math import *
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages 


############################
# ERROR ON THE ABSOLUTE MAGNITUDE VERSUS REDSHIFT
def errorZ(fl,mess,ax):

   # step en mag
   zstep=0.1
   zstart=0
   znbstep=70

   medpVec=[]
   mednVec=[]
   zVec=[]
   
   #set the axis
   ax.axis([0,6,0,4])

   # Loop over the redshift bins
   for rm in xrange(znbstep) :

      #limit inf and sup
      zlimi=zstart+zstep*rm
      zlims=zstart+zstep*(rm+1)

      #new condition with the magnitude range 
      condA = (zp>zlimi) & (zp<=zlims) & (eval("filt"+str(fl))>0) & (eval("eabsmag"+str(fl))>-0.01) & (eval("eabsmag"+str(fl))<4)
 
      # be sure to have galaxies in a sufficient number
      if len( zp[condA]) > 3 :
       
       #compute the median value
       medp=np.median(eval("eabsmag"+str(fl))[condA])
       
       # Add the result to the vector
       medpVec.append(medp)
       zVec.append(zstart+zstep*(rm+0.5))
      
   #plot the error versus z
   ax.plot(zVec,medpVec, color='b',linestyle='-')

   # labels
   ax.annotate(mess,xy=(zstart+0.2,3.5),color="black", fontsize=15)


############################
# Compare rest-frame colors 
def colcol(zmin,zmax,flt1,flt2,flt3,mess,ax):

    #set the axis
   ax.axis([-1,7,-1,7])

   # Color computed from the absolute magnitude
   colAbsmag1=eval("absmag"+str(flt1) + " -  absmag"+str(flt2))
   colAbsmag2=eval("absmag"+str(flt2) + " -  absmag"+str(flt3))

   #plot color versus color
   condA = (zp>zmin) & (zp<=zmax)
   ax.scatter(col1[condA],colAbsmag1[condA], s=1, color='b',alpha=0.5,marker='s')
   ax.scatter(col2[condA],colAbsmag2[condA], s=1, color='r',alpha=0.5,marker='s')
 
   #droite
   xy=np.array([-10,10])
   ax.plot(xy,xy)


############################
# Compare rest-frame colors 
def william(zlimi,flt1,flt2,flt3,mess,ax):

    #set the axis
   ax.axis([-2,3,-1,7])

   # Color computed from the absolute magnitude
   colAbsmagX=eval("absmag"+str(flt2) + " -  absmag"+str(flt3))
   colAbsmagY=eval("absmag"+str(flt1) + " -  absmag"+str(flt2))
 
   #plot color versus color
   condA = (zp>zlimi) & (zp<=zlimi+1)
   ax.scatter(colAbsmagX[condA],colAbsmagY[condA], s=1, color='b',alpha=0.5,marker='s')

   ax.scatter(col2[condA],col1[condA], s=1, color='r',alpha=0.5,marker='s')



############################
# Compare absolute magnitude old versus new
def absOld(zlimi,flt1,flt2,flt3,mess,ax):

    #set the axis
   ax.axis([-27,-13,-27,-13])

   #droite
   ax.plot(xy,xy)

   # condition
   condA = (zp>zlimi) & (zp<=zlimi+1)

   # old versus new absolute magnitudes
   absmagX=eval("absmag"+str(flt1))
   absmagY=eval("absmagOld"+str(flt1))
   ax.scatter(absmagX[condA],absmagY[condA], s=1, color='b',alpha=0.5,marker='s')

   absmagX=eval("absmag"+str(flt2))
   absmagY=eval("absmagOld"+str(flt2))
   ax.scatter(absmagX[condA],absmagY[condA], s=1, color='r',alpha=0.5,marker='s')

   absmagX=eval("absmag"+str(flt3))
   absmagY=eval("absmagOld"+str(flt3))
   ax.scatter(absmagX[condA],absmagY[condA], s=1, color='g',alpha=0.5,marker='s')


############################
# Compare chi2 old versus new
def chiOldNew(zlimi,mess,ax):

    #set the axis
   ax.axis([0,1000,0,1000])

   #droite
   ax.plot(xy,xy)

   # condition
   condA = (zp>zlimi) & (zp<=zlimi+1)

   # old versus new absolute magnitudes
   ax.scatter(chiOld[condA],chi[condA], s=1, color='b',alpha=0.5,marker='s')




############################
# FILTER USED FOR THE ABSOLUTE MAGNITUDE VERSUS REDSHIFT
def fltZ(fl,mess,ax):
   
   #set the axis
   ax.axis([-2,30,0,5.])
   flhist=np.arange(300)*0.25-3

   #new condition with the magnitude range 
   condA = (zp>2) & (zp<=3)
   # be sure to have galaxies in a sufficient number
   if len( zp[condA]) > 3 :              
     #plot the error versus z
     ax.hist((eval("filt"+str(fl))[condA])-0.25, bins=flhist, histtype='stepfilled', normed=True, color='b', label=mess)

   #new condition with the magnitude range 
   condA = (zp>3) & (zp<=4)
   # be sure to have galaxies in a sufficient number
   if len( zp[condA]) > 3 :              
     #plot the error versus z
     ax.hist((eval("filt"+str(fl))[condA]), bins=flhist, histtype='stepfilled', normed=True, color='g', label=mess)

   #new condition with the magnitude range 
   condA = (zp>4) & (zp<=5)
   # be sure to have galaxies in a sufficient number
   if len( zp[condA]) > 3 :              
     #plot the error versus z
     ax.hist((eval("filt"+str(fl))[condA])+0.25, bins=flhist, histtype='stepfilled', normed=True, color='r')
     # labels
     ax.annotate(mess,xy=(4,4),color="black", fontsize=15)



########
# READ THE INPUT FILE

# Read the first argument with the name of the Le Phare output catalogue
fileIn='magabs.out'
fileInOld='../fortran/magabs.out'
catIn=open(fileIn,'r')
catInOld=open(fileInOld,'r')
print "Name of the Le Phare output catalogue : ",fileIn
print "Name of the Le Phare output catalogue Old : ",fileInOld


# Loop over the filters
nbFilt=39
mst=""
idmst=""
mstOld=""
idmstOld=""
# create the string to read the apparent mag
for i in xrange(nbFilt) :
    mst=mst+",mag"+str(i)
    idmst=idmst+","+str(i+24)
    mstOld=mstOld+",magOld"+str(i)
    idmstOld=idmstOld+","+str(i+24)
magst=""
idmagst=""
magstOld=""
idmagstOld=""
# create the string to read the absolute mag
for i in xrange(nbFilt) :
    magst=magst+",absmag"+str(i)
    idmagst=idmagst+","+str(i+141)
    magstOld=magstOld+",absmagOld"+str(i)
    idmagstOld=idmagstOld+","+str(i+141)
# create the string to read the absolute mag errors
emagst=""
idemagst=""
for i in xrange(nbFilt) :
    emagst=emagst+",eabsmag"+str(i)
    idemagst=idemagst+","+str(i+180)
# create the string to read the filter used
filt=""
idfilt=""
filtOld=""
idfiltOld=""
for i in xrange(nbFilt) :
    filt=filt+",filt"+str(i)
    idfilt=idfilt+","+str(i+219)
    filtOld=filtOld+",filtOld"+str(i)
    idfiltOld=idfiltOld+","+str(i+180)


# Extract from the ascii file
commandst = "Id,zp,chi"+mst+magst+emagst+filt+",ssfr_inf,ssfr_med,ssfr_sup,col1_inf,col1,col1_sup,col2_inf,col2,col2_sup= np.loadtxt(catIn, dtype='float', usecols=(0,1,11"+idmst+idmagst+idemagst+idfilt+",279,280,281,282,283,284,285,286,287), unpack=True )"
commandstOld = "IdOld,zpOld,chiOld"+mstOld+magstOld+filtOld+"= np.loadtxt(catInOld, dtype='float', usecols=(0,1,11"+idmstOld+idmagstOld+idfiltOld+"), unpack=True )"

print commandst
print commandstOld

# transform the string into a command
exec(commandst)
exec(commandstOld)


###############################
######### ERRORS VERSUS Z


# All the figures will be collected in a single pdf file 
pdfOut = PdfPages('figures.pdf')


## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$z$', ha='center')
f.text(0.04, 0.5, '$absolute \; magnitude \; uncertainties$', va='center', rotation='vertical')

# Define the subplots and pass the panel
errorZ(31,"NUV",axarr[0,0])  
errorZ(3,"R",axarr[0,1])   
errorZ(8,"J",axarr[1,0])  
errorZ(10,"K",axarr[1,1]) 

# store the figure in a PDF
savefig(pdfOut,format='pdf')

## Create a figure with an array with nbRowZ*nbColZ subpanels
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$absolute \; magnitude \; filters$', ha='center')
f.text(0.04, 0.5, '$N$', va='center', rotation='vertical')

# Define the subplots and pass the panel
fltZ(31,"NUV",axarr[0,0])  
fltZ(3,"R",axarr[0,1])   
fltZ(8,"J",axarr[1,0])  
fltZ(10,"K",axarr[1,1]) 

# store the figure in a PDF
savefig(pdfOut,format='pdf')


## Create a figure with an array with 2x2 subpanels color-color plot
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$color new$', ha='center')
f.text(0.04, 0.5, '$color old$', va='center', rotation='vertical')

# Define the subplots and pass the panel
colcol(0,1,31,3,10,"0<z<1",axarr[0,0])  
colcol(1,2,31,3,10,"1<z<2",axarr[0,1])  
colcol(2,3,31,3,10,"2<z<3",axarr[1,0])  
colcol(0,4,31,3,10,"0<z<4",axarr[1,1])  


# store the figure in a PDF
savefig(pdfOut,format='pdf')



## Create a figure with an array with 2x2 subpanels color-color plot
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$R-K$', ha='center')
f.text(0.04, 0.5, '$NUV-R$', va='center', rotation='vertical')

# Define the subplots and pass the panel
william(0,31,3,10,"0<z<1",axarr[0,0])  
william(1,31,3,10,"1<z<2",axarr[0,1])  
william(2,31,3,10,"2<z<3",axarr[1,0])  
william(3,31,3,10,"3<z<4",axarr[1,1])  


##########
# OLD verus NEW mag abs
#########

## Create a figure with an array with 2x2 subpanels color-color plot
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$M old$', ha='center')
f.text(0.04, 0.5, '$M new$', va='center', rotation='vertical')
   
xy=np.array([-100,10])

# Define the subplots and pass the panel
absOld(0,31,3,10,"0<z<1",axarr[0,0])  
absOld(1,31,3,10,"1<z<2",axarr[0,1])  
absOld(2,31,3,10,"2<z<3",axarr[1,0])  
absOld(3,31,3,10,"3<z<4",axarr[1,1])  


# store the figure in a PDF
savefig(pdfOut,format='pdf')

##########
# OLD verus NEW chi2
#########

## Create a figure with an array with 2x2 subpanels color-color plot
clf()
f, axarr = subplots(2,2,sharex=True,sharey=True)
# No space between the figures
subplots_adjust( hspace=0,wspace=0)

# label of the figure
f.text(0.5, 0.04, '$chi2 old$', ha='center')
f.text(0.04, 0.5, '$chi2 new$', va='center', rotation='vertical')
   
xy=np.array([-100,10000])

# Define the subplots and pass the panel
chiOldNew(0,"0<z<1",axarr[0,0])  
chiOldNew(1,"1<z<2",axarr[0,1])  
chiOldNew(2,"2<z<3",axarr[1,0])  
chiOldNew(3,"3<z<4",axarr[1,1])  


# store the figure in a PDF
savefig(pdfOut,format='pdf')

pdfOut.close()

print 'end'


