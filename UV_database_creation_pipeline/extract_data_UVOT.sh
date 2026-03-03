#!/usr/bin/bash
j=0
datapath="../../../Data"

echo "Extracting the required columns for UVOT Catalogue"
stilts tpipe ifmt=fits \
cmd='addcol -after DEC POSERR sqrt(square(RA_ERR*cosDeg(DEC))+square(DEC_ERR))' \
cmd='keepcols "OBSID SRCNUM RA DEC POSERR \
UVW2_FLUX UVW2_FLUX_ERR UVM2_FLUX UVM2_FLUX_ERR \
UVW1_FLUX UVW1_FLUX_ERR U_FLUX U_FLUX_ERR \
B_FLUX B_FLUX_ERR V_FLUX V_FLUX_ERR UVW2_QUALITY_FLAG \
UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
V_QUALITY_FLAG UVW2_EXTENDED UVM2_EXTENDED UVW1_EXTENDED \
U_EXTENDED B_EXTENDED V_EXTENDED"' \
omode=out ofmt=fits in=$datapath/Swift-UVOT/uvotssc1_1.fit#1 \
out=$datapath/Swift-UVOT/temp.fits

echo "Dividing the time table in various part"
python Utilities/divide_time_table.py "Swift-UVOT" "uvotssc1_1.fit"

echo "Match of each time table with the main catalogue"
j=0
for f in $( ls $datapath/Swift-UVOT/temp_part*.fits); do
    stilts tmatch2 ifmt1=fits ifmt2=fits in1="$datapath/Swift-UVOT/temp.fits" \
    in2=$f matcher=exact \
    values1="OBSID" values2="OBSID" \
    ocmd="addcol MJD_START isoToMjd(DATE_MIN)" \
    ocmd='keepcols "SRCNUM RA DEC POSERR \
    UVW2_FLUX UVW2_FLUX_ERR UVM2_FLUX UVM2_FLUX_ERR \
    UVW1_FLUX UVW1_FLUX_ERR U_FLUX U_FLUX_ERR \
    B_FLUX B_FLUX_ERR V_FLUX V_FLUX_ERR UVW2_QUALITY_FLAG \
    UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
    V_QUALITY_FLAG UVW2_EXTENDED UVM2_EXTENDED UVW1_EXTENDED \
    U_EXTENDED B_EXTENDED V_EXTENDED MJD_START"' \
    join=1and2 find=best1 ofmt=csv out=$datapath/Swift-UVOT/temp_match.csv \
    progress=none

    if [ $j -gt 0 ]; then
  		echo "Concatenating matched slices"
	    stilts tcat ifmt=csv in="$datapath/Swift-UVOT/UVOT_complete.csv \
        $datapath/Swift-UVOT/temp_match.csv" \
        icmd="replacecol SRCNUM toInteger(SRCNUM)" \
        out=$datapath/Swift-UVOT/temp.csv
		    
	    mv $datapath/Swift-UVOT/temp.csv $datapath/Swift-UVOT/UVOT_complete.csv
    else
	    mv $datapath/Swift-UVOT/temp_match.csv $datapath/Swift-UVOT/UVOT_complete.csv
    fi
    let j=j+1
done

echo "Generating the RA DEC file"
stilts tpipe ifmt=csv cmd='keepcols "SRCNUM RA DEC POSERR"' \
cmd='sort SRCNUM' omode=out ofmt=fits in=$datapath/Swift-UVOT/UVOT_complete.csv \
out=$datapath/Swift-UVOT/uvot_ra_dec.fits

echo "Generating the Entries file"
stilts tpipe ifmt=csv cmd='keepcols "SRCNUM UVW2_FLUX UVW2_FLUX_ERR UVM2_FLUX UVM2_FLUX_ERR \
UVW1_FLUX UVW1_FLUX_ERR U_FLUX U_FLUX_ERR \
B_FLUX B_FLUX_ERR V_FLUX V_FLUX_ERR UVW2_QUALITY_FLAG \
UVM2_QUALITY_FLAG UVW1_QUALITY_FLAG U_QUALITY_FLAG B_QUALITY_FLAG \
V_QUALITY_FLAG UVW2_EXTENDED UVM2_EXTENDED UVW1_EXTENDED \
U_EXTENDED B_EXTENDED V_EXTENDED MJD_START"' \
cmd='sort SRCNUM' \
omode=out ofmt=csv in=$datapath/Swift-UVOT/UVOT_complete.csv \
out=$datapath/Swift-UVOT/Entries_UVOT_csv

rm $datapath/Swift-UVOT/temp.fits $datapath/Swift-UVOT/temp_*.fits $datapath/Swift-UVOT/temp_match.csv $datapath/Swift-UVOT/UVOT_complete.csv



