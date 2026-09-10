#!/usr/bin/bash

# Path to the output data and the Galex slices
# Note that the existence of a sub-repository named Galex
datapath="../../../Data"
datapath_ext="/media/julien/TOSHIBA EXT/Copie_non_miroir/Master/TFM/Data/Galex slices"
echo "Random match count" > $datapath/Galex/count_random_match.txt

j=0
# Loop over the Galex galactic slices
for f in "$datapath_ext"/GUVCat_*.csv.gz; do

    # Extract the required columns from the galex slices
	echo "Extracting the required columns for " "$f"
	stilts tpipe ifmt=csv cmd='keepcols "objid ra dec fuv_flux \
	fuv_fluxerr nuv_flux nuv_fluxerr nuv_poserr"' \
	omode=out ofmt=csv in="$f" out=$datapath/Galex/temp.csv

    # Copy the temporary file for the all possible matches
    cp $datapath/Galex/temp.csv $datapath/Galex/temp_all.csv

    # Count the number of matches with the UVOT-OM catalogue
	echo "Performing match with UVOT-OM catalogue for " "$f"
	COUNT="$(stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
	in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	values1="ra dec maxReal(3.6,3*nuv_poserr)" values2="RA DEC 3*POSERR" \
	join=1and2 find=best params=1 omode=count progress=none \
    | cut -d ":" -f 3)"

    echo "Number of matches : ${COUNT}"

    # If the number of matches is greater than 0
    if [ "${COUNT}" -gt 0 ]; then

        # Perform the match with the best point creating a new angDist column and selecting 
        # the kept columns
        # Note that the match considers the max of 3*sigma and 3 times the mean sigma over 
        # the whole distribution for Galex point
	    stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
	    in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	    values1="ra dec maxReal(3.6,3*nuv_poserr)" values2="RA DEC 3*POSERR" \
        ocmd="addcol angDist skyDistanceDegrees(ra_1,dec_1,RA_2,DEC_2)*3600" \
        ocmd='keepcols "objid ra_1 dec_1 fuv_flux fuv_fluxerr nuv_flux \
	    nuv_fluxerr nuv_poserr SRCNUM_CUV SRCNUM_OM SRCNUM_UVOT RA_2 DEC_2 POSERR COUNT angDist"' \
	    join=1and2 find=best ofmt=csv params=1 out=$datapath/Galex/temp_match.csv \
        progress=none

        # Perform the match with all points creating a new angDist column and selecting 
        # the kept columns
        # Note that the match considers the max of 3*sigma and 3 times the mean sigma over 
        # the whole distribution for Galex point
	    stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp_all.csv \
	    in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	    values1="ra dec maxReal(3.6,3*nuv_poserr)" values2="RA DEC 3*POSERR" \
        ocmd="addcol angDist skyDistanceDegrees(ra_1,dec_1,RA_2,DEC_2)*3600" \
        ocmd='keepcols "objid ra_1 dec_1 fuv_flux fuv_fluxerr nuv_flux \
	    nuv_fluxerr nuv_poserr SRCNUM_CUV SRCNUM_OM SRCNUM_UVOT RA_2 DEC_2 POSERR COUNT angDist"' \
	    join=1and2 find=all ofmt=csv params=1 out=$datapath/Galex/temp_match_all.csv \
        progress=none

        # Compute the match with random values moved from 1 arcmin
        echo "Computing random match values"
        stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
        in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
        values1="ra dec maxReal(3.6,3*nuv_poserr)" values2="RA+0.01667 DEC 3*POSERR" \
        join=1and2 find=best params=1 omode=count \
        progress=none >> $datapath/Galex/count_random_match.txt

        # If the slice is not the first one
	    if [ $j -gt 0 ]; then

            # Concatenate the results for the best matches
      		echo "Concatenating matched slices"
		    stilts tcat ifmt=csv in="$datapath/Galex/match_Galex_UVOT_OM_best_max36.csv \
            $datapath/Galex/temp_match.csv" \
            icmd="replacecol SRCNUM_CUV toInteger(SRCNUM_CUV)" \
            icmd="replacecol SRCNUM_OM toInteger(SRCNUM_OM)" \
            icmd="replacecol SRCNUM_UVOT toInteger(SRCNUM_UVOT)" \
            out=$datapath/Galex/temp.csv

            # Concatenate the results for all matches
		    stilts tcat ifmt=csv in="$datapath/Galex/match_Galex_UVOT_OM_all_max36.csv \
            $datapath/Galex/temp_match_all.csv" \
            icmd="replacecol SRCNUM_CUV toInteger(SRCNUM_CUV)" \
            icmd="replacecol SRCNUM_OM toInteger(SRCNUM_OM)" \
            icmd="replacecol SRCNUM_UVOT toInteger(SRCNUM_UVOT)"\
            out=$datapath/Galex/temp_all.csv

            # Move the temporary files to the concatenated results file
		    mv $datapath/Galex/temp.csv $datapath/Galex/match_Galex_UVOT_OM_best_max36.csv
            mv $datapath/Galex/temp_all.csv $datapath/Galex/match_Galex_UVOT_OM_all_max36.csv
            
	    # Otherwise
        else
            # Move the temporary files to the concatenated results file directly
		    mv $datapath/Galex/temp_match.csv $datapath/Galex/match_Galex_UVOT_OM_best_max36.csv
            mv $datapath/Galex/temp_match_all.csv $datapath/Galex/match_Galex_UVOT_OM_all_max36.csv
	    fi
	    let j=j+1

    else
        echo "No match for " $f
    fi
    
    echo ""
done
