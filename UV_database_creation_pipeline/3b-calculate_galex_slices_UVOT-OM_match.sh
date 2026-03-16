#!/usr/bin/bash
j=0
datapath="../../../Data"
echo "Random match count" > $datapath/Galex/count_random_match.txt
for f in $( ls $datapath/Galex/GUVCat_*.csv.gz); do

	echo "Extracting the required columns for " $f
	stilts tpipe ifmt=csv cmd='keepcols "objid ra dec fuv_flux \
	fuv_fluxerr nuv_flux nuv_fluxerr nuv_poserr"' \
	omode=out ofmt=csv in=$f out=$datapath/Galex/temp.csv

    cp $datapath/Galex/temp.csv $datapath/Galex/temp_all.csv

	echo "Performing match with UVOT-OM catalogue for " $f
	COUNT="$(stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
	in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	values1="ra dec 3*nuv_poserr" values2="RA DEC 3*POSERR" \
	join=1and2 find=best params=1 omode=count progress=none \
    | cut -d ":" -f 3)"

    echo "Number of matches : ${COUNT}"

    if [ "${COUNT}" -gt 0 ]; then
	    stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
	    in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	    values1="ra dec 3*nuv_poserr" values2="RA DEC 3*POSERR" \
        ocmd="addcol angDist skyDistanceDegrees(ra_1,dec_1,RA_2,DEC_2)*3600" \
        ocmd='keepcols "objid ra_1 dec_1 fuv_flux fuv_fluxerr nuv_flux \
	    nuv_fluxerr nuv_poserr SRCNUM_OM SRCNUM_UVOT RA_2 DEC_2 POSERR COUNT angDist"' \
	    join=1and2 find=best ofmt=csv params=1 out=$datapath/Galex/temp_match.csv \
        progress=none

	    stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp_all.csv \
	    in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
	    values1="ra dec 3*nuv_poserr" values2="RA DEC 3*POSERR" \
        ocmd="addcol angDist skyDistanceDegrees(ra_1,dec_1,RA_2,DEC_2)*3600" \
        ocmd='keepcols "objid ra_1 dec_1 fuv_flux fuv_fluxerr nuv_flux \
	    nuv_fluxerr nuv_poserr SRCNUM_OM SRCNUM_UVOT RA_2 DEC_2 POSERR COUNT angDist"' \
	    join=1and2 find=all ofmt=csv params=1 out=$datapath/Galex/temp_match_all.csv \
        progress=none

        echo "Computing random match values"
        stilts tmatch2 ifmt1=csv ifmt2=csv in1=$datapath/Galex/temp.csv \
        in2=$datapath/SOURCE_complete_galaxies.csv matcher=skyerr \
        values1="ra dec 0" values2="RA+0.01667 DEC 3*POSERR" \
        join=1and2 find=best params=1 omode=count \
        progress=none >> $datapath/Galex/count_random_match.txt

	    if [ $j -gt 0 ]; then
      		echo "Concatenating matched slices"
		    stilts tcat ifmt=csv in="$datapath/Galex/match_Galex_UVOT_OM_best.csv \
            $datapath/Galex/temp_match.csv" \
            icmd="replacecol SRCNUM_OM toInteger(SRCNUM_OM)" \
            icmd="replacecol SRCNUM_UVOT toInteger(SRCNUM_UVOT)" \
            out=$datapath/Galex/temp.csv

		    stilts tcat ifmt=csv in="$datapath/Galex/match_Galex_UVOT_OM_all.csv \
            $datapath/Galex/temp_match_all.csv" \
            icmd="replacecol SRCNUM_OM toInteger(SRCNUM_OM)" \
            icmd="replacecol SRCNUM_UVOT toInteger(SRCNUM_UVOT)"\
            out=$datapath/Galex/temp_all.csv
		    
		    mv $datapath/Galex/temp.csv $datapath/Galex/match_Galex_UVOT_OM_best.csv
            mv $datapath/Galex/temp_all.csv $datapath/Galex/match_Galex_UVOT_OM_all.csv
	    else
		    mv $datapath/Galex/temp_match.csv $datapath/Galex/match_Galex_UVOT_OM_best.csv
            mv $datapath/Galex/temp_match_all.csv $datapath/Galex/match_Galex_UVOT_OM_all.csv
	    fi
	    let j=j+1

    else
        echo "No match for " $f
    fi
    
    echo ""
done
