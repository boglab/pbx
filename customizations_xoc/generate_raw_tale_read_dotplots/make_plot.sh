fileprefix=`basename "${1}" | cut -c 11- | rev | cut -c 7- | rev`
/opt/MUMmer3.23/nucmer --maxmatch --nosimplify -l 10 -c 10 --maxgap 150 -p ${fileprefix} ${1} ${1}
/opt/MUMmer3.23/delta-filter -l 50 ${fileprefix}.delta > ${fileprefix}-filtered.delta
/opt/MUMmer3.23/mummerplot -t postscript -p ${fileprefix} ${fileprefix}-filtered.delta
