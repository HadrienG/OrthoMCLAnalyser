#!/bin/bash

# requirements: mcl, orthomcl, blast mysql
# this program configures an runs the orthomcl program.

Usage="$(basename "$0") [-h] [-n database Name] [-u database login name] [-p database password] [-g Genome list] [-d .gb Directory] [-o Output directory] -- RunOrthoMCL

    -h show this useful help
    -n <str> Name of the MySQL database
    -u <str> Your database login name
    -p <str> Your database password
    -g <txt> GenomeList.txt. See the example file for syntax
    -d <dir> Directory containing your genbank files
    -o <dir> Output directory"

while getopts ":hn:u:p:g:d:o:" option
do
  case $option in
    h) echo "$Usage"
    exit 0
      ;;
    n) DBName=$OPTARG
      ;;
    u) DBUser=$OPTARG
      ;;
    p) DBPass=$OPTARG
      ;;
    g) GenomeList=$OPTARG
      ;;
    d) GenomeDir=$OPTARG
      ;;
    o) OutDir=$OPTARG
    if [[ ! -d $OutDir ]]
      then
        echo "Creating Output Directory"
        mkdir -p $Outdir
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

# Check for requirements
command -v mysql >/dev/null 2>&1 || { echo >&2 "MySQL is not installed or is not in the PATH. Aborting\n"; exit 1; }
command -v orthomclInstallSchema >/dev/null 2>&1 || { echo >&2 "orthoMCL is not installed or is not in the PATH. Aborting\n"; exit 1; }
command -v blastp >/dev/null 2>&1 || { echo >&2 "blastp is not installed or is not in the PATH. Aborting.\n"; exit 1; }
command -v mcl >/dev/null 2>&1 || { echo >&2 "mcl is not installed or is not in the PATH. Aborting\n"; exit 1; }
command -v perl >/dev/null 2>&1 || { echo >&2 "perl is not installed or is not in the PATH. Aborting\n"; exit 1; }
command -v python >/dev/null 2>&1 || { echo >&2 "python is not installed or is not in the PATH. Aborting\n"; exit 1; }
perl -MDBI -e ";" >/dev/null 2>&1 || { echo >&2 "the DBI libraries are not installed. Aborting\n"; exit 1; }
python -c "from Bio import SeqIO" >/dev/null 2>&1 || { echo >&2 "Biopython is not installed. Aborting\n"; exit 1; }


# Check if the parameters are correctly set
if [ -z ${DBName+x} ]; then echo "$Usage" >&2; printf "\nPlease provide a database name\n" >&2; exit 1; else echo "Database Name is set to '$DBName'"; fi
if [ -z ${DBUser+x} ]; then echo "$Usage" >&2; printf "\nPlease provide a database login name\n" >&2; exit 1; else echo "Database Login Name is set to '$DBUser'"; fi
if [ -z ${DBPass+x} ]; then echo "$Usage" >&2; printf "\nPlease provide a database password\n" >&2; exit 1; else echo "Database Login Name is set to '$DBPass'"; fi
if [ -z ${GenomeList+x} ]; then echo "$Usage" >&2; printf "\nPlease provide a genome list\n" >&2; exit 1; else echo "Genome List is set to '$GenomeList'"; fi
if [ -z ${GenomeDir+x} ]; then echo "$Usage" >&2; printf "\nPlease provide a input directory\n" >&2; exit 1; else echo "Genome Directory is set to '$GenomeList'"; fi
if [ -z ${Outdir+x} ]; then echo "$Usage" >&2; printf "\nPlease provide an output directory\n" >&2; exit 1; else echo "Output Directory is set to '$Outdir'"; fi

# Create orthomcl.config
cat >$OutDir/orthomcl.config << EOL
# orthomcl.config generated automatically by RunOrthoMCL.sh
dbVendor=mysql
dbConnectString=dbi:mysql:$DBName:mysql_local_infile=1
dbLogin=$DBUser
dbPassword=$DBPass
similarSequencesTable=SimilarSequences
orthologTable=Ortholog
inParalogTable=InParalog
coOrthologTable=CoOrtholog
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
EOL
echo "$OutDir/orthomcl.config created"

# Install the database. Usage = orthomclInstallSchema config_file sql_log_file table_suffix
mysql -u $DBUser -p$DBPass -e "create database if not exists $DBName";
orthmclInstallSchema $OutDir/orthomcl.config $OutDir/InstallSchema.log
echo "Database $DBName and schema created"

# Conversion: gb to faa. Usage = python gbtofaa.py -i input -o output.
# orthomclAdjustFasta. Usage = orthomclAdjustFasta taxon_code fasta_file id_field
mkdir -p $OutDir/compliantFasta
while IFS=$' ' read genome code field
  do
    python Python/gbtofaa.py -i $GenomeDir$genome -o $OutDir/${genome%.*}
    cd $OutDir/compliantFasta
    orthomclAdjustFasta $code $OutDir/${genome%.*}.faa $field
    cd $WorkDir
  done <${GenomeList}
echo "Compliant fasta files created"