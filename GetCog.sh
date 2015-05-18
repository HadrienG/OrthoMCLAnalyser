#!/bin/bash

# requirements: rpsblast, biopython
# this program blast the .faa of the core and pan Directories against the
# Cog database

Usage="$(basename "$0") [-h] [-i Input Directory] [-o Output Directory] -- GetCog.sh

    -h show this useful help
    -i <dir> Input Directory
    -o <dir> Output directory"

while getopts ":hi:o:" option
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
if [ -z "${Input+x}" ]; then echo "$Usage" >&2; printf "\nThe Input Directory should be your RunOrthoMCL output.\n" >&2; exit 1; else echo "Input Directory is set to '$Input'"; fi
if [ -z "${OutDir+x}" ]; then echo "$Usage" >&2; printf "\nPlease provide an output directory\n" >&2; exit 1; else echo "Output Directory is set to '$OutDir'"; fi


# Download the CoG database
if [ ! -f $OutDir/Database/Cog.rps ]
  then
    cd $OutDir
    mkdir -p Database && cd Database
    curl -O ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
    tar xzf Cog_LE.tar.gz
    cd ..
fi
DbName="Database/Cog"
cd $WorkDir
if [ ! -d $OutDir/CogInfo ]
  then
  cd $OutDir
  mkdir -p CogInfo && cd CogInfo
  curl -O ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab
  curl -O ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab
  cd ..
fi

# Blast the Core and Pan genomes against the Cog db
cd $WorkDir/$OutDir
mkdir -p $WorkDir/$OutDir/Xml
for fasta in ../$Input/Infos/Genomes/*.faa
do
  rpsblast -query "$fasta" -db $DbName -evalue 0.00001 -outfmt 5 -out Xml/$(basename $fasta).xml
  if [ "$?" != 0 ]
    then
      echo "rpsblast failed. Aborting." >&2
      exit 1
  fi
done


# Parse the blast results
cd $WorkDir/$OutDir/Xml
python "$ScriptPath"/Python/parserpsblast.py
cd $WorkDir/$OutDir
python "$ScriptPath"/Python/cog2plot.py
cd $WorkDir

# Get the ratio of COG hits / total number of proteins
TotalCore=$(grep -ch ">" $Input/Infos/Genomes/*Core.faa | awk '{ SUM += $1 } END { print SUM }')
TotalPan=$(grep -ch ">" $Input/Infos/Genomes/*Pan.faa | awk '{ SUM += $1 } END { print SUM }')
COGCore=$(wc -l $OutDir/CoreCogList.txt | awk '{ print $1 }')
COGPan=$(wc -l $OutDir/PanCogList.txt | awk '{ print $1 }')
UnknownCore=$((TotalCore-COGCore))
UnknownPan=$((TotalPan-PanCore))

# Produces Plots of the functional annotations
cd $WorkDir
Rscript "$ScriptPath"/R/MeanCoreVsPan.R -i $Input/InfosSummary.txt -o $OutDir
Rscript "$ScriptPath"/R/Categories.R -i $OutDir -o $OutDir -c $UnknownCore -p $UnknownPan
