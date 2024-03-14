#----------------------------------------------------------
#  ZPHOT
#    Makefile
#----------------------------------------------------------
#
#  UNIT
#
#       File:	Makefile 
#
#       Author:	Benoit Vandame    :  Many Thanks guy !!
#       Updated: 19/12/2013
#     
#----------------------------------------------------------
#
#  DESCRIPTION
#
#       Makefile for LE PHARE  Programs
#
#----------------------------------------------------------
#
#
################################################
#######       Compilation      #################
################################################
#
#
OS := $(shell uname -s)
ifeq ($(OS), Darwin)
    # For macOS, brew install gcc
    FC         =   g++-13
    CFLAGS     = -g -std=c++17 -fopenmp -Wall -fno-inline -Wno-psabi
endif
ifeq ($(OS), Linux)
    # For Linux
    FC         = g++
    #CFLAGS     = -g -std=c++11 -O4 -fopenmp -Wall -fno-inline
    CFLAGS     = -ggdb -std=c++17 -Ofast -fopenmp -Wall -fno-inline -march=native -mtune=native -mno-vzeroupper -funroll-loops -floop-unroll-and-jam -fomit-frame-pointer -finline -mavx2 -ftree-vectorize  #-fopt-info-vec-all
endif
LIBS       = ##/usr/local/lib/libMinuit2.a
#
########       Create the LIBRARY      ##########
LIB_LEPHARE  = lephare_lib.a
AR         = ar
RANLIB     = ranlib
#
########         Source archive        ##########
ARCHIVE = lephare++_v0.0.tar
#
########        Compiling source       ##########

FILES  = \
        filter\
        filter_extinc\
        sedtolib\
        mag_gal\
        zphota

OBJECTS = \
        globals.o\
        keyword.o\
        flt.o\
        SED.o\
	SEDLib.o\
        ext.o\
        opa.o\
        mag.o\
        oneElLambda.o\
        cosmology.o\
        onesource.o\
        PDF.o\
	photoz_lib.o\

all: $(FILES)
	if [ ! -d "../../bin" ]; then mkdir ../../bin; fi
	cp $(FILES) ../../bin/

%.o: %.cpp *h 
	$(FC) $(CFLAGS) -c  $*.cpp  -o $*.o   

%: %.cpp $(LIB_LEPHARE) *h
	$(FC) $(CFLAGS) $*.cpp -o $*  $(LIB_LEPHARE)  $(LIBS)


$(FILES): $(LIB_LEPHARE)

$(LIB_LEPHARE): $(OBJECTS)
	rm -f $(LIB_LEPHARE)
	$(AR) rv $(LIB_LEPHARE) $(OBJECTS)
	$(RANLIB) $(LIB_LEPHARE)
#
#
#################################################
########            Cleaning             ########
#################################################
#
#
clean:
	rm -f \#*
	rm -f *.o
	rm -f *.a
	rm -f *~

	rm -f minuit/\#*
	rm -f minuit/*.o
	rm -f minuit/*.a
	rm -f minuit/*~
	rm -f ../../bin/*

#
#
#################################################
########            Archiving           #########
#################################################
#
#
archi:
	rm -f $(ARCHIVE)
	tar cvf $(ARCHIVE) *.cpp *.h
	gzip $(ARCHIVE)
	echo "Archive: $(ARCHIVE).gz"
#
#
################################################
#######     Create working directory        ####
################################################
work:
	mkdir -p $(LEPHAREWORK)
	mkdir -p $(LEPHAREWORK)/filt
	mkdir -p $(LEPHAREWORK)/lib_bin
	mkdir -p $(LEPHAREWORK)/lib_mag

doc:
	doxygen ../doc/Doxyfile

