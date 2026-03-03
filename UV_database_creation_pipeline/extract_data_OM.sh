#!/usr/bin/bash
j=0
datapath="../../../Data"

echo "Extracting the required columns for OM Catalogue"
stilts tpipe ifmt=fits cmd='keepcols "OBSID SRCNUM RA DEC POSERR \
UVW2_AB_FLUX UVW2_AB_FLUX_ERR UVM2_AB_FLUX UVM2_AB_FLUX_ERR \
UVW1_AB_FLUX UVW1_AB_FLUX_ERR U_AB_FLUX U_AB_FLUX_ERR \
B_AB_FLUX B_AB_FLUX_ERR V_AB_FLUX V_AB_FLUX_ERR UVW2_QUALITY_FLAG \
UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
V_QUALITY_FLAG UVW2_EXTENDED_FLAG UVM2_EXTENDED_FLAG UVW1_EXTENDED_FLAG \
U_EXTENDED_FLAG B_EXTENDED_FLAG V_EXTENDED_FLAG"' \
omode=out ofmt=fits in=$datapath/XMM-Newton/XMM-OM-SUSS6.2.fits#1 \
out=$datapath/XMM-Newton/temp.fits

echo "Dividing the time table in various part"
python Utilities/divide_time_table.py "XMM-Newton" "XMM-OM-SUSS6.2.fits"

echo "Match of each time table with the main catalogue"
j=0
for f in $( ls $datapath/XMM-Newton/temp_part*.fits); do
    stilts tmatch2 ifmt1=fits ifmt2=fits in1="$datapath/XMM-Newton/temp.fits" \
    in2=$f matcher=exact \
    values1="OBSID" values2="OBSID" \
    ocmd='keepcols "SRCNUM RA DEC POSERR \
    UVW2_AB_FLUX UVW2_AB_FLUX_ERR UVM2_AB_FLUX UVM2_AB_FLUX_ERR \
    UVW1_AB_FLUX UVW1_AB_FLUX_ERR U_AB_FLUX U_AB_FLUX_ERR \
    B_AB_FLUX B_AB_FLUX_ERR V_AB_FLUX V_AB_FLUX_ERR UVW2_QUALITY_FLAG \
    UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
    V_QUALITY_FLAG UVW2_EXTENDED_FLAG UVM2_EXTENDED_FLAG UVW1_EXTENDED_FLAG \
    U_EXTENDED_FLAG B_EXTENDED_FLAG V_EXTENDED_FLAG MJD_START"' \
    join=1and2 find=best1 ofmt=csv out=$datapath/XMM-Newton/temp_match.csv \
    progress=none

    if [ $j -gt 0 ]; then
  		echo "Concatenating matched slices"
	    stilts tcat ifmt=csv in="$datapath/XMM-Newton/OM_complete.csv \
        $datapath/XMM-Newton/temp_match.csv" \
        icmd="replacecol SRCNUM toInteger(SRCNUM)" \
        out=$datapath/XMM-Newton/temp.csv
		    
	    mv $datapath/XMM-Newton/temp.csv $datapath/XMM-Newton/OM_complete.csv
    else
	    mv $datapath/XMM-Newton/temp_match.csv $datapath/XMM-Newton/OM_complete.csv
    fi
    let j=j+1
done

echo "Generating the RA DEC file"
stilts tpipe ifmt=csv cmd='keepcols "SRCNUM RA DEC POSERR"' \
cmd='sort SRCNUM' omode=out ofmt=fits in=$datapath/XMM-Newton/OM_complete.csv \
out=$datapath/XMM-Newton/SUSS6_ra_dec.fits

echo "Generating the Entries file"
stilts tpipe ifmt=csv cmd='keepcols "SRCNUM UVW2_AB_FLUX \
UVW2_AB_FLUX_ERR UVM2_AB_FLUX UVM2_AB_FLUX_ERR \
UVW1_AB_FLUX UVW1_AB_FLUX_ERR U_AB_FLUX U_AB_FLUX_ERR \
B_AB_FLUX B_AB_FLUX_ERR V_AB_FLUX V_AB_FLUX_ERR UVW2_QUALITY_FLAG \
UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
V_QUALITY_FLAG UVW2_EXTENDED_FLAG UVM2_EXTENDED_FLAG UVW1_EXTENDED_FLAG \
U_EXTENDED_FLAG B_EXTENDED_FLAG V_EXTENDED_FLAG MJD_START"' \
cmd='sort SRCNUM' \
omode=out ofmt=csv in=$datapath/XMM-Newton/OM_complete.csv \
out=$datapath/XMM-Newton/Entries_OM_csv

rm $datapath/XMM-Newton/temp.fits $datapath/XMM-Newton/temp_*.fits $datapath/XMM-Newton/temp_match.csv $datapath/XMM-Newton/OM_complete.csv



