#!/bin/bash

# requirements: rpsblast, biopython
# this program blast the .faa of the core and pan Directories against the
# Cog database

Usage="$(basename "$0") [-h] [-i Input Directory] [-o Output Directory] -- GetCog.sh

    -h show this useful help
    -i <dir> Input Directory
    -o <dir> Output directory"

while getopts ":ho:" option
do
  case $option in
    h) echo "$Usage"
    exit 0
      ;;
    i) Input=$OPTARG
      ;;
    o) OutDir=$OPTARG
    if [[ ! -d $OutDir ]]
      then
        echo "Creating Output Directory"
        mkdir -p "$OutDir"
    fi
      ;;
    \?) printf "illegal option -%s\n" "$OPTARG" >&2
    echo "$Usage" >&2
    exit 1
      ;;
  esac
done
shift $((OPTIND-1))
WorkDir=$(pwd)
pushd "$(dirname "$0")" >/dev/null
ScriptPath=$(pwd)
popd > /dev/null


# Check for requirements
command -v python >/dev/null 2>&1 || { printf "python is not installed or is not in the PATH. Aborting\n" >&2; exit 1; }
python -c "from Bio import SeqIO" >/dev/null 2>&1 || { printf "Biopython is not installed. Aborting\n" >&2; exit 1; }
command -v rpsblast >/dev/null 2>&1 || { printf "rpsblast is not installed or is not in the PATH. Aborting\n" >&2; exit 1; }


# Check if the parameters are correctly set
if [ -z "${OutDir+x}" ]; then echo "$Usage" >&2; printf "\nPlease provide an output directory\n" >&2; exit 1; else echo "Output Directory is set to '$OutDir'"; fi


# Download the CoG database
if [ ! -f /$OutDir/Database/Cog.rps ]
  then
    cd $OutDir
    mkdir -p Database
    curl -O ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
    tar xzf Cog_LE.tar.gz
fi
DbName="$OutDir/Database/Cog"
if [ ! -d /$OutDir/CogInfo ]
  then
  cd $OutDir
  mkdir -p CogInfo
  curl -O ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab
  curl -O ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab


# Blast the Core and Pan genomes against the Cog db
mkdir -p $OutDir/Xml
for fasta in $Input/Infos/Genomes/*.faa
do
  rpsblast -i "$fasta" -d .$DbName -e 0.00001 -m 7 -a 8 -o $OutDir/Xml/"$fasta".xml
done


# Parse the blast results
cd $OutDir/Xml
python "$ScriptPath"/Python/parserpsblast.py
cd $OutDir
python "$ScriptPath"/Python/cog2plot.py


# Produces Plots of the functional annotations
cd $OutDir
Rscript "$ScriptPath"/R/MeanCoreVsPan.R -i $Input/InfosSummary.txt -o $OutDir
Rscript "$ScriptPath"/R/Categories.R -i InputDir -o $OutDir









exit 0
