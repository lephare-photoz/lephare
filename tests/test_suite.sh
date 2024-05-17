#!/usr/bin/env bash

export OMP_NUM_THREADS='30'
tmpdir=$(python -c "import lephare as lp; print(f'{lp.LEPHAREDIR}')" | tail -n 1) # use tail to ignore the on-import print statement
tmpwork=$(python -c "import lephare as lp; print(f'{lp.dm.LEPHAREWORK}')" | tail -n 1)
export LEPHAREDIR=$tmpdir
export LEPHAREWORK=$tmpwork

echo "LEPHAREDIR:" $LEPHAREDIR # TODO remove these at the end
echo "LEPHAREWORK:" $LEPHAREWORK

# TODO check the paths in the files themselves - as in, the list files and such

# Get the data
#curl -s -o ${LEPHAREDIR}/COSMOS.in https://raw.githubusercontent.com/lephare-photoz/lephare-data/main/examples/COSMOS.in
echo "STARTING DOWNLOAD"
if [ -d "$LEPHAREDIR" ]; then
    curl -o ${LEPHAREDIR}/examples/COSMOS.in https://raw.githubusercontent.com/lephare-photoz/lephare-data/main/examples/COSMOS.in
    if [ $? -eq 0 ]; then
        echo "Download successful."
    else
        echo "Download failed."
    fi
else
    echo "Directory $LEPHAREDIR does not exist."
fi

# Variables
export CONFIG_FILE="$LEPHAREDIR/examples/COSMOS.para"

# Test suite commands
filter -c $CONFIG_FILE

# sedtolib -c ./COSMOS.para -t S --STAR_SED $LEPHAREDIR/examples/STAR_MOD_ALL.list --LIB_ASCII YES
# mag_gal -c ./COSMOS.para -t S --LIB_ASCII YES --STAR_LIB_OUT ALLSTAR_COSMOS

# sedtolib -c $LEPHAREDIR/examples/COSMOS.para -t Q --QSO_SED  $LEPHAREDIR/sed/QSO/SALVATO09/AGN_MOD.list
# mag_gal -c $LEPHAREDIR/examples/COSMOS.para -t Q --MOD_EXTINC 0,1000  --EB_V 0.,0.1,0.2,0.3 --EXTINC_LAW SB_calzetti.dat --LIB_ASCII NO  --Z_STEP 0.04,0,6 --LIB_ASCII YES

# sedtolib -c $LEPHAREDIR/examples/COSMOS.para -t G --GAL_SED $LEPHAREDIR/examples/COSMOS_MOD.list  --GAL_LIB LIB_VISTA
# mag_gal  -c $LEPHAREDIR/examples/COSMOS.para -t G --GAL_LIB_IN LIB_VISTA --GAL_LIB_OUT VISTA_COSMOS_FREE --MOD_EXTINC 18,26,26,33,26,33,26,33  --EXTINC_LAW SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat  --EM_LINES EMP_UV  --EM_DISPERSION 0.5,0.75,1.,1.5,2. --Z_STEP 0.04,0,6 --LIB_ASCII YES

# cat_out=zphot_short.out

# zphota -c $LEPHAREDIR/examples/COSMOS.para --CAT_IN $LEPHAREDIR/examples/COSMOS.in --CAT_OUT $cat_out --ZPHOTLIB VISTA_COSMOS_FREE,ALLSTAR_COSMOS,QSO_COSMOS  --ADD_EMLINES 0,100 --AUTO_ADAPT YES   --Z_STEP 0.04,0,6 --CAT_LINES 1,100 --SPEC_OUT YES --PARA_OUT $LEPHAREDIR/examples/output.para --VERBOSE NO --ZFIX NO --PDZ_OUT $LEPHAREWORK/zphota/

# python $LEPHAREDIR/examples/figuresLPZ.py $cat_out
# python $LEPHAREDIR/examples/spec.py *.spec -d pdf -o $LEPHAREWORK/zphota/spec

# mv $cat_out Id*.spec figuresLPZ.pdf $LEPHAREWORK/zphota/