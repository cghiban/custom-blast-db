#!/bin/sh

#set -x

### run it in parallel like this:
#
# parallel --jobs 3 bin/phy_blast_get_accn_data.sh {} ::: fasta-20180306/200-[01]*.txt
#
###

INFILE=$1
FORMAT=$2
if [[ -z $INFILE ]]; then
    echo "E: missing argument <file>"
    exit 1
fi

if [[ ! -f $INFILE ]]; then
    echo "E: $INFILE not found!"
    exit 1
fi

if [[ -z $FORMAT ]]; then
    FORMAT=fasta
fi


echo "#-------------------------------"
echo -e " + processing $INFILE"

OUTFILE=${INFILE/\.*/.$FORMAT}
#if [[ ! $OUTFILE =~ \.fasta$ ]]; then
#    OUTFILE="$OUTFILE.fasta"
#fi

if [[ -f $OUTFILE ]]; then
    #echo " - output file ${OUTFILE} already exists. Skipping..."
    exit 0
fi

echo " + output to be stored in $OUTFILE"

IDLIST=$(cat $INFILE | grep -v '^$' |perl -pe 's/\n/,/'|perl -pe 's/,+$//')

if [[ -z $IDLIST ]]; then
    echo " - not IDs given. continuing..."
    exit 0
fi
#echo " + IDLIST=${IDLIST}"
#exit 0

EPOSTOUT=$(curl -sk -X POST -d "id=${IDLIST}&email=dnalcadmin@cshl.edu&db=nuccore" https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi)

#set -x

#echo " + EPOSTOUT=${EPOSTOUT}"
#echo -e " + EPOSTOUT=${EPOSTOUT}"

WE=$(echo -e $EPOSTOUT | grep WebEnv | perl -pwe 's#.*<WebEnv>(.*?)</WebEnv>.*#$1#g')

#echo " + WE=${WE}"

if [[ -z $WE ]]; then
    echo -e "Can't extract WebEnv from EPOSTOUT:\n${EPOSTOUT}"
    exit 1
fi

# get the data in fasta/gb format

curl -s -o $OUTFILE "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=${FORMAT}&WebEnv=${WE}&query_key=1&retstart=0&retmax=500"

CEXIT=$?
if [[ $CEXIT != 0 ]]; then
    echo " E: curl exit code: ${CEXIT}"
fi

if [[ $FORMAT == "fasta" ]]; then
    NUMSEQ=$(grep "^>" $OUTFILE|wc -l)
    echo -e " = ${OUTFILE}\t${NUMSEQ}"
else
    ENTRIES=$(grep "^LOCUS" $OUTFILE|wc -l)
    echo -e " = ${OUTFILE}\t${ENTRIES}"
fi

sleep 1
