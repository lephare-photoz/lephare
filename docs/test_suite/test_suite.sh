#!/usr/bin/env bash
# This script gives an example of using the command line arguments 
# We call some Python code for establishing the directories and
# downloading data required for the example but if you have your own files 
# present they would not be required.

# Set the number of threads for OpenMP if your system requires
# export OMP_NUM_THREADS='30'

# Set the LEPHAREDIR and LEPHAREWORK variables automatically to defaults
export LEPHAREDIR=$(python -c "import lephare as lp; print(f'{lp.LEPHAREDIR}')" | tail -n 1) # use tail to ignore the on-import print statement
export LEPHAREWORK=$(python -c "import lephare as lp; print(f'{lp.dm.LEPHAREWORK}')" | tail -n 1) # use tail to ignore the on-import print statement
# Alternatively set them to your preferred local locations
# export LEPHAREDIR=/my/local/lepharedir
# export LEPHAREWORK=/my/local/lepharework

# Get the required data files and example scripts for the run using the Python functionality
# All files downloaded from https://github.com/lephare-photoz/lephare-data/
python -c "import lephare as lp; \
lp.data_retrieval.get_auxiliary_data( \
keymap=lp.default_cosmos_config, \
additional_files=[ \
'examples/COSMOS.in', \
'examples/COSMOS.para', \
'examples/STAR_MOD_ALL.list', \
'sed/QSO/SALVATO09/AGN_MOD.list', \
'examples/COSMOS_MOD.list', \
'examples/figuresLPZ.py', \
'examples/spec.py', \
])"

# Set catalog and config file names
CAT_IN="$LEPHAREDIR/examples/COSMOS.in"
CAT_OUT="zphot_short.out"
CONFIG_FILE="$LEPHAREDIR/examples/COSMOS.para"

# Run commands
echo "Running filter..."
filter -c $CONFIG_FILE

echo "Running sedtolib and mag_gal for stars..."
sedtolib -c $CONFIG_FILE \
    -t S \
    --STAR_SED $LEPHAREDIR/examples/STAR_MOD_ALL.list \
    --LIB_ASCII YES
mag_gal -c $CONFIG_FILE \
    -t S \
    --LIB_ASCII YES \
    --STAR_LIB_OUT ALLSTAR_COSMOS

echo "Running sedtolib and mag_gal for QSOs..."
sedtolib -c $CONFIG_FILE \
    -t Q \
    --QSO_SED $LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list
mag_gal -c $CONFIG_FILE \
    -t Q \
    --MOD_EXTINC 0,1000  \
    --EB_V 0.,0.1,0.2,0.3 \
    --EXTINC_LAW SB_calzetti.dat \
    --LIB_ASCII NO \
    --Z_STEP 0.04,0,6 \
    --LIB_ASCII YES

echo "Running sedtolib and mag_gal for galaxies..."
sedtolib -c $CONFIG_FILE \
    -t G \
    --GAL_SED $LEPHAREDIR/examples/COSMOS_MOD.list  \
    --GAL_LIB LIB_VISTA 
mag_gal -c $CONFIG_FILE \
    -t G \
    --GAL_LIB_IN LIB_VISTA \
    --GAL_LIB_OUT VISTA_COSMOS_FREE \
    --MOD_EXTINC 18,26,26,33,26,33,26,33  \
    --EXTINC_LAW SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat  \
    --EM_LINES EMP_UV \
    --EM_DISPERSION 0.5,0.75,1.,1.5,2. \
    --Z_STEP 0.04,0,6 \
    --LIB_ASCII YES

echo "Running zphota..."
zphota -c $CONFIG_FILE \
    --CAT_IN $CAT_IN \
    --CAT_OUT $CAT_OUT \
    --ZPHOTLIB VISTA_COSMOS_FREE,ALLSTAR_COSMOS,QSO_COSMOS  \
    --ADD_EMLINES 0,100 \
    --AUTO_ADAPT YES \
    --Z_STEP 0.04,0,6 \
    --CAT_LINES 1,100 \
    --SPEC_OUT YES \
    --PARA_OUT $LEPHAREDIR/examples/output.para \
    --VERBOSE NO \
    --ZFIX NO \
    --PDZ_OUT $LEPHAREWORK/zphota/

echo "Generating figures and specs..."
python $LEPHAREDIR/examples/figuresLPZ.py $CAT_OUT
python $LEPHAREDIR/examples/spec.py *.spec -d pdf -o $LEPHAREWORK/zphota/spec

echo "Moving output files..."
mv $CAT_OUT Id*.spec figuresLPZ.pdf $LEPHAREWORK/zphota/
