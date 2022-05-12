#!/usr/bin/env bash

## Copyright (c) 2020 Translational Genomics Research Institute
##
## This software may be modified and distributed under the terms
## of the MIT license.  See the LICENSE file for details.
##
## Major Contributors: Christophe Legendre
## Minor Contributors: Jessica Aldrich

### @@@@@@@@
### COMMENTS
### tgen.mutation_burden.sh
### number of cpus minimum required: 5 due to the piping and -@ in samtools view
### this tool is specific for a run with Tumor_versus_Constitutional use case, with somatic WES or WGS
### Inputs required:
### 1) a BED file with Targets regions
### 2) a VCF file with somatic calls between the tumor and constitutional Bam files
### 3) a constitutional Bam file
### 4) a tumor bam file
###
### Default values are seen in the Default Values section in the function called getOptions below but the default can be modified using appropriate options
###
### flags for samtools explained: 0x10 is << read reverse strand >> ; 0x400 is << read is PCR or optical duplicate >>
### @@@@@@@@



##@@@@@@@@@@@@@
## FUNCTIONS
##@@@@@@@@@@@@@
#function clean_old_files(){
#	echo -e "cleaning old files if exist ..."
#	for F in $*
#	do
#		if [[ ${F} == "" || ${F} == " " ]]; then continue ; fi
#		echo -e "deleting file ... ${F}" ;
#		rm --verbose -f "${F}"
#	done
#	sleep 1
#}

function pipelinize(){
	## output style requested by internal pipeline developer for JSON import
	local LOUTFILE="${1}"
	local PIPELINE="${2}"
	if [[ ${PIPELINE} == "yes" ]]
	then
		cp "${LOUTFILE}" ${LOUTFILE}.temp
		echo -e "making output ala pipeline-style for json import"
		( awk 'NR<8' ${LOUTFILE}.temp ; awk 'NR>=8' ${LOUTFILE}.temp | awk '{FS=OFS="\t" ; print $1,$2,$3,$13,$14,$15}' ) > ${LOUTFILE}
		if [[ $? == 0 ]] ; then rm ${LOUTFILE}.temp ; fi
	fi
}

function add_file_to_list(){
	local L=$1 ; local F=$2
	if [[ ${#L[@]} -ne 0  ]] ; then L="${L[@]} ${F}" ; else L="${F}" ; fi
	echo "${L[@]}"
}

function check_Rscript(){
	type Rscript >/dev/null 2>&1 || { echo >&2 "Require  \"Rscript\"   but it's not in the PATH.  Aborting."; exit 1; }
}

function check_samtools_version(){
	local min_samtools_version_required=1.7
	type samtools >/dev/null 2>&1 || { echo >&2 "Require  \"samtools\"   but it's not in the PATH.  Aborting."; exit 1; }
}

function check_bgzip_tabix(){
	type bgzip >/dev/null 2>&1 || { echo >&2 "Require  \"bgzip\"   but it's not in the PATH.  Aborting."; exit 1; }
	type tabix >/dev/null 2>&1 || { echo >&2 "Require  \"tabix\"   but it's not in the PATH.  Aborting."; exit 1; }
}


function check_bedtools(){
	type bedtools >/dev/null 2>&1 || { echo >&2 "Require  \"bedtools\"   but it's not in the PATH.  Aborting."; exit 1; }
}

function usage(){
	echo -e "##\n##\n## OPTIONS [M == Mandatory, O = Optional ]"
	echo -e "##
<---- mandatory options ---->
--bed \t\tBED file representing the TARGETS Regions [M]
--vcf \t\tVCF file with somatic variants [M]
--nbam \t\tConstitutional or Normal BAM file [M]
--tbam \t\tCase or Tumor BAM file [M]
\n<---- optional ----->
--mindepth \tMinimum read depth for both BAMS at each positions --> Callable Space [O, Default value is 10]
--min-read-per-strand \tMinimum read depth for each Strand for Case/Tumor BAM at each positions [O, Default value is 1]
--min-base-qual \tMinimum base quality to keep the read in samtools depth cmd [O, Default value is 5]
--min-map-qual \tMinimum Mapping quality to keep the read in samtools view cmd [O, Default value is 5]
--outprefix \tProvide a prefix name to default output filenames [O, default is NULL]
--pqarms \tprovide a 2-columns FILE with column 1 as the coordinates (1:0-121535434) of the arms and the column 2 as the name of the arms (such as p or q); see file examples; WARNING: if you have lots of small regions in that file, the plot will be very long to make and may be not easy to read or not well formatted; The plot has been optimized for regions similar to regions of p and q arms chromosomes in size; up to 5 regions per chromosome might be the maximum for good plot read
--keep-unflt-outfile \tIf you need to keep the unfiltered samtools depth's output file, use this option  [O, Default is 'no' ]
--skip=plot \tThis option is only used when --pqarms is provided and the user does not want to get a plot with ALL the regions given; Note; For each region given, a bar will exist in plot; so if you have thousands of region in BED file, you will get thousands of bar in the plot which may take a while to make; If used, then you may consider this as running in the regions in parallel without getting the plot at the end
--samtools-depth-only \tAdvanced_users_Only; It defies the purpose of this tool as it skips the calculation of the 'Mutation Burden' and only runs the 'samtools depth' commands
--force \tIf the contigs between the two bam files do not perfectly match both at the name and size levels, the script will stop; use this option to force the run knowing that might generate wrong calculi. Use at your own risk unless you are sure that the contigs difference will not impact the Mutation Burden calculation, such as the difference may be in alternate contigs not present in the VCF; [O, Default is 'no']
\n<---- advanced users options  ----->
--region \tYou may provide a region formatted as follow 1:100000-200000 ; coordinates value are inclusive here; this will speed up the output while testing if you modified the current script ; the region will only be used with the samtools view to limit the searching region space in bams ; the BED file is still mandatory when using this option ; also help for DEBUGGING by reducing the space;  [O, Default is NULL]
--region-name \t  Only used if --pqarms is provided ; For Internal usage when recursive call made [O, Default is NULL]
--do-not-add-dbsnp-pct-to-plot \t enable plotting the percentage percentage of dbsnp in given regions ; required --pqarms option with a 3-column file provided to --pqarms option\t  [O, Default is \"no\"]
--flag-dbsnp \t Default value hardcoded is \"RS\" which means that we count the RSID present in column three (3) of the VCF; if vcf is annotated with another flag for knowing if variant is present or not in DBSNP, please provide the FLAG This Flag should be defined in INFO column; (if not ask for updating the script by creating an issue in github); Example --flag-dbsnp DBSNP152 [O, default is RS ]
-j|--gnu-parallel-jobs \t This is for more advance users; this option will require LOTS of CPUS available if value if greater than 1 ; see README to learn about how many cpus total will be require [O, Default is 1 ]
-t|--threads \tfor DEBUGGING; Use with CAUTION ; If you want to save 5 minutes to your run you may use upto 2 threads per samtools command ; this will therefore increase the need of cpus to 10 cpus to run this script; Do not Use unless necessary ; [O, Default is 1]
--verbose \tFor DEBUGGING purposes; print to sdtout the content of the stats outfile [O, Default is 'no']
\n<---- internal option ---->
--pct-dbsnp  Only used if --pqarms is provided and if file from option --pqarms contains a third column with dbsnp percentage -- AND -- this option is only used internally  \t  [O, Default is NULL]
\n<---- help ----->
-h|--help \tprint the options for usage and exit
"

echo -e "Example:
## 1) normal usage
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM}
## 2) normal usage with extra options
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --prefix 'anyStringHere_withOnlyAlphaNumericCharacter'
## 3) if you want to keep the output from samtools depth comamnd, use --keep-unflt-outfile option
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --keep-unflt-outfile option
## 4) if you want to use different min depth values
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --keep-unflt-outfile option --mindepth 20 --min-read-per-strand 2 --min-base-qual 20 --min-map-qual 20
## 5) if you want to use only a specific region without provinding a bed file
$0 --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --region 1:1000000-2000000
## 6) if you want to get the mutation Burden per p and q arms of the chromosomes, provide a 2-column tabulated file with column #1 as arm's coordinate and column #2 as the name of that region
$0 --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --pqarms \${MY_PQ_ARM_FILE}
## 7) same as Example #6 with dbsnp percentage not added to plot
$0 --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --pqarms \${MY_PQ_ARM_FILE} --do-not-add-dbsnp-pct-to-plot
## 8) debugging usage, add verbose option to your command
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --region 1:1000000-2000000 --verbose
## 9) for TGen pipeline, additional options are needed
$0 --bed \${MYBED} --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --sample MYSAMPLENAME --library KKKKKK --rg MyRG --pipeline --verbose --outfile TMB_SAMPLE_STATS.txt

## ** Advanced Users **
## ** --threads and --gnu-parallel usage **
## Warning: a run with threads=1 and gnu-parallel-jobs=1 will ALREADY consume 4 cpus at 100% ; so calculate consequences; for instance: 4x1x1=4 ; 4x2x2=16 ; 4x1x4=16
## Here with example #9, example of running `basename $0` in an optimize way with trade off between speed and runtime [ ~ 35 minutes ]
## 9) same as Example #7 with gnu-parallel-jobs enabled
$0 --vcf \${MYVCF} --nbam \${MYNBAM} --tbam \${MYTBAM} --pqarms \${MY_PQ_ARM_FILE} --do-not-add-dbsnp-pct-to-plot --threads 1 --gnu-parallel-jobs 4
"
}


function checkEV(){
	local ev=$1
	local msg="$2"
	if [[ ${ev} -ne 0 ]]
	then
		echo -e "\nERROR: ${msg} FAILED; Aborting!\n" ;
		exit 1 ;
	fi
}


function checkIfInt(){
	local I=$1
	if [[ ! ${I} =~ ^[0-9]*$ ]] ; then echo -e "\nERROR: Positive Integer Number between 1 and infinity is Expected;
	Found:  << ${I} >> ; Aborting!\n" ; exit 1 ; fi
}

function checkFile(){
	local F=$1
	if [[ ! -e ${F} || ! -f ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >>; Aborting!" ; exit 1 ; fi ;
}

function checkDir(){
	local D=$1
	if [[ ! -e ${D} || ! -d ${D} ]] ; then echo -e "DIR NOT FOUND << ${D} >>; Aborting!" ; exit 1 ; fi ;
}

function init_defaults(){
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## Default Values
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	MINDEPTH=10						## Apply to the SUM of Reads from both strand
	MINREADSPERSTRAND=1		## Specifically Applied to Positions in Case/Tumor BAM
	MIN_BASE_QUAL=5
	MIN_MAP_QUAL=5
	THREADS=1
	GNU_PARALLEL=1
	FORCE="no"
	KEEP_OUTFILE_UNFILTCOUNT_FSRS="no"
	RUN_ONLY_DEPTH_CAPTURE="no"
	VERBOSE="no"
	PREFIX=""
	PADDED_BED=""
	PQARMS_FILE=""
	PQARM_COORD_NAME=""
	REGION_NAME="NA"
	FLAG_DBSNP="RS"
	PCT_DBSNP="NA"
	SAMPLE=""
	LIBRARY=""
	READ_GROUP=""
	SKIP_PLOT="no"
	OUTFILE=""
	MAX_CALLABLE_SPACE=0
	PIPELINE="no"
	DIR_TEMP=""
	ADD_DBSNP_TO_PLOT="yes" ; 		## it was no by default b/c bug found while calculating percentage; yes for test
	REGION_INTERVAL=""  					## FOR DEBUGGING or Tests purposes ONLY: we can use a region such as  10:102798692-102801121
	LIST_ERROR_FILES="err.f0x10_F0x400.Constitutional.$$.log err.F0x10_F0x400.Constitutional.$$.log err.f0x10_F0x400.Case.$$.log err.F0x10_F0x400.Case.$$.log"
	for F in ${LIST_ERROR_FILES[@]} ; do rm "${F}" &>/dev/null ; done
}

function getOptions(){

	## init the defaults values before capturing the user values that may override the defaults based on the used options
	init_defaults

	if ! options=`getopt -o hj:t: -l bed:,vcf:,nbam:,tbam:,mindepth:,min-read-per-strand:,min-base-qual:,min-map-qual:,outprefix:,sample:,library:,rg:,outfile:,force,keep-unflt-outfile,verbose,threads:,gnu-parallel-jobs:,region:,region-name:,pqarms:,pct-dbsnp:,dir-temp:,flag-dbsnp:,do-not-add-dbsnp-pct-to-plot,samtools-depth-only,skip-plot,pipeline,help -- "$@"`
	then
	# something went wrong, getopt will put out an error message for us
		echo -e "<--- ERROR in Arguments ----> \n" ; usage ; exit 1 ;
	fi
	eval set -- "$options"
	while [[ $# -gt 0 ]]
	do
		# echo -e "$1 === $2"  ; ## Uncomment this lien for debugging
		# getopt recognize the params and put unrecognized params after '--', so '--' is a very important flag for latter work.
		# for options with a required argument, an additional shift is required
		case $1 in
		--bed) PADDED_BED=$2 ; checkFile ${2};  shift ;;
		--vcf) VCF=${2} ; checkFile ${2}; shift ;;
		--nbam) NBAM=${2}  ; checkFile ${2}; shift ;;
		--tbam) TBAM=${2}  ; checkFile ${2}; shift ;;
		--pqarms) PQARMS_FILE=${2} ; checkFile ${2} ; shift ;;
		--mindepth) MINDEPTH=${2}  ; checkIfInt ${2};  shift ;;
		--min-read-per-strand) MINREADSPERSTRAND=${2} ; checkIfInt ${2};  shift ;;
		--min-base-qual) MIN_BASE_QUAL=${2}  ; checkIfInt ${2};  shift ;;
		--min-map-qual) MIN_MAP_QUAL=${2}  ; checkIfInt ${2};  shift ;;
		--outprefix) PREFIX=${2} ; shift ;;
		--region) REGION_INTERVAL=${2} ; shift ;;
		--flag-dbsnp) FLAG_DBSNP="$2" ; shift ;;
		--pct-dbsnp) PCT_DBSNP=${2} ; shift ;;
		--do-not-add-dbsnp-pct-to-plot) ADD_DBSNP_TO_PLOT="no" ;;
		--region-name) REGION_NAME=$( echo ${2} | sed 's/ \+/_/g' ) ; shift ;;
		--keep-unflt-outfile) KEEP_OUTFILE_UNFILTCOUNT_FSRS="yes" ;;
		--samtools-depth-only) RUN_ONLY_DEPTH_CAPTURE="yes" ; KEEP_OUTFILE_UNFILTCOUNT_FSRS="yes" ;;
		--force) FORCE="yes" ;;
		--verbose) VERBOSE="yes" ;;
		--sample) SAMPLE="$2" ; shift ;;
		--library) LIBRARY="$2" ; shift ;;
		--rg) READ_GROUP="$2" ; shift ;;
		--skip-plot) SKIP_PLOT="yes"  ;;
		--outfile) OUTFILE="$2" ; shift ;;
		--pipeline) PIPELINE="yes" ;; ## special option request by pipeline developer BT ;-)
		--dir-temp) DIR_TEMP="$2" ; shift ;;
		-t|--threads) THREADS=$2 ; checkIfInt ${2} ; shift ;;
		-j|--gnu-parallel-jobs) GNU_PARALLEL=${2} ; checkIfInt ${2} ; shift ;;
		-h|--help) usage ; exit ;;
		(--) shift ; echo "--" ;;
		(-*) break ; echo -e "$0: error - unrecognized option $1\n\n" 1>&2 ; usage;  exit 1  ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ; usage;  exit 1  ;;
		esac
		shift
	done

	## we check the PREFIX before using it for the output filename
	if [[ ${PREFIX} != "" ]] ; then
		if [[ ${PREFIX} =~ ^[a-zA-Z0-9\_\-\.]*$ && ! ${PREFIX} =~ [?@\#\$\&\(\){}]+ ]] ; then
		    if [[ ! ${PREFIX} =~ [_\.]$ ]] ; then PREFIX="${PREFIX}." ; fi ;
		else
		    echo -e "\nERROR: Unwanted Character in the Prefix String provided; Please only use AlphaNumerical character and absolutely NO special ones\ngiven prefix was: ${PREFIX}" ; exit 1 ;
		fi
	fi
	## we check if the index file is present
	echo "VCF == ${VCF} in dir currdir: ${PWD}"
	tabix --list-chrom ${VCF} &>/dev/null
	if [[ $? -ne 0 ]]; then echo -e "ERROR: Missing VCF index; please use command 'bcftools index --tbi ${VCF}' ; Aborting" ; exit 1 ; fi
	## we check if at least ONE of BED or REGION is provided
	if [[ ${PADDED_BED} == "" && ${REGION_INTERVAL} == "" ]] ; then echo -e "\nERROR: You need to provide at least one of the two options --bed or --region; Found NONE; Aborting." ; exit 1 ; fi
	if [[ ${PQARMS_FILE} != "" ]] ; then check_Rscript ; fi
	## we check if threads greater than 3; due to the implementation, having more than 2 threads in useless due to the limitation of samtools depth
	 if [[ ${THREADS} -gt 4 ]] ; then THREADS=4 ; fi ## HARDCODED VALUE ; if someone demonstrates to me that more than 4 cpus per samtools view speeds thing up, I will change that HARDCODED VALUE of 4 to any other suggested value;
	##capture basename for adding filename to header file
	NBAM_BN=$(basename ${NBAM}) ;
	TBAM_BN=$(basename ${TBAM} ) ;
	## Defining Output filenames here based on Case BAM name and PREFIX if given
	OUTFILE_STATS=${PREFIX}${TBAM_BN}.coverage.mutation_burden.stats.txt ;
	OUTFILE_UNFILTCOUNT_FSRS=${PREFIX}${TBAM_BN}.CallablePositions.unflt.txt.gz
	if [[ -e ${OUTFILE_STATS} && ${FORCE} == "no" ]] ; then echo -e "WARNING: file << ${OUTFILE_STATS} >> already EXISTS; Please Modify your prefix or delete that file first, or use --force to overwrite the file; Aborting;" ; exit 1 ; elif [[ -e ${OUTFILE_STATS} && ${FORCE} == "yes"  ]] ; then echo -e "WARNING: as you have used --force option, file << ${OUTFILE_STATS} >> which if EXISTED has been OVERWRITTEN." ; fi
	## print ALL inputs and Variables
	LIST_VARS=( PADDED_BED VCF NBAM TBAM MINDEPTH MINREADSPERSTRAND MIN_BASE_QUAL MIN_MAP_QUAL PREFIX OUTFILE_STATS OUTFILE_UNFILTCOUNT_FSRS THREADS KEEP_OUTFILE_UNFILTCOUNT_FSRS VERBOSE REGION_INTERVAL PQARMS_FILE )
	if [[ "${VERBOSE}" == "yes" ]] ; then for ITEM in "${LIST_VARS[@]}" ; do echo -e "${ITEM}\t==\t${!ITEM}" ; done ; echo "" ; fi
	ORIGINAL_COMMAND="$(echo -e "$(realpath $0) \c") $( (($#)) && printf ' %q' "$@")"
}



############
### MAIN ###
############
echo $(date)
## checking if executable in PATH
check_samtools_version
check_bgzip_tabix
check_bedtools

## capturing user given options to current script
getOptions $@

## for capturing the intermediate errors in piped commands
set -eu -o pipefail


## @@@@@@@@@@@@@@@@@
## PREPROCESSING
## @@@@@@@@@@@@@@@@@
echo -e "WORKING_DIRECTORY is: ${PWD} "

## check bam or cram integrity ; this will minimize the errors with the next samtools view commands
samtools quickcheck -qv ${NBAM} ${TBAM} > bad_bams.$$.fofn && ( echo "--> ALL BAMS OK <-- " ; rm bad_bams.$$.fofn) || ( echo -e "ERROR: Some files failed check, see file bad_bams.fofn" ; exit 1 )

## check if Bam Header @SQ are identical ; i.e. Checking if same contigs in all bams [both at name and size levels] (JK's request)
if [[ $(paste <( samtools view -H ${NBAM} | grep -E "^@SQ" | cut -f1-3 | sort | sed 's/[ \+\t]/__/g' ) <(samtools view -H ${TBAM} | grep -E "^@SQ" | cut -f1-3 | sort | sed 's/[ \+\t]/__/g' ) | tee contigs.${NBAM_BN}_${TBAM_BN}.$$.txt | awk -F'\t' '{ if($1!=$2) { print $0} } ' | wc -l ) -ne 0 ]] ;
then
	echo -e "Contigs are different from one file to the other ; check your input BAM files; The BAM files you have been processed the same way using the same reference genome. see file contigs.${NBAM_BN}_${TBAM_BN}.$$.txt " ;
	if [[ ${FORCE} == "no" ]] ; then exit 1 ; else echo -e "WARNING: as you have used --force option, you decided to continue even though the reference genome used was apparently different since the contigs differ; this may lead to unexpected results about the mutation Burden." ; fi
else
	echo -e "--> SAME CONTIGS OK  <-- " ; rm contigs.${NBAM_BN}_${TBAM_BN}.$$.txt  ;
fi



## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## By COORDINATE REGIONS
## We check it the user provided a file with coordinates for the p and q arms [NOTE: actually it could be any file with a list coordinates; 1 coordinate per line]
## All the coordinates in this file MUST also be represented in the BED Region Input File; If not, wrong results will arise.
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if [[ ${PQARMS_FILE} != "" && ${REGION_INTERVAL} == "" ]] ;
then
	## TODO: Rethink the entire approach for parallel; this old style is ugly but it works
	## init some variables for the current "if" section
	LOSF="" ## List OUT STATS FILE
	LFAILURE="" ## list of Failure
	PREFIX_PQARM=""
	LOG_PQARM_FAILURE="log.${PREFIX}${TBAM_BN}.coverage.mutation_burden.stats.txt.list_pqarms_failures.txt"

	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## PARALLEL
	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if [[ ${GNU_PARALLEL} -lt 2 ]]
	then
		## SERIAL PROCESSING since GNU_PARALLEL == 1 or 0
		## we process the file with coordinates that should be in column one with this format CHRNAME:START-END
		while read line
		do
			LOPT="" ## ListOPTion for tgen_mutation_burden script
			if [[ ${PADDED_BED} != "" ]] ; then LOPT="${LOPT[@]} --bed ${PADDED_BED}" ; fi
			if [[ ${FORCE} == "yes" ]] ; then LOPT="${LOPT[@]} --force " ; fi

			COORD="$(echo "${line}" | cut -f1 )" ## this column is mandatory
			PQARM_COORD_NAME="$(echo "${line}" | cut -f2 )" ## this column is mandatory

			## PCT_DBSNP="$(echo "${line}" | cut -f3 )" ## this column is optional ; actually it can be any covariate you want, not necessary dbsnp pct but anything else related to that \$COORD
			PREFIX_PQARM=${PREFIX}pqarms.$(echo "${COORD}" | sed 's/[:-]/_/g')
			OSF=${PREFIX_PQARM}.${TBAM_BN}.coverage.mutation_burden.stats.txt ;  ##OUTPUT STATS FILE  ; There is an Issue here b/c I do not see why it needs a dot to work. ## TO BE TESTED MORE
			if [[ "${VERBOSE}" == "yes" ]] ; then echo -e "VERBOSE: ${COORD} ${PQARM_COORD_NAME} ${PREFIX_PQARM} ${OSF}" ; fi

			## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			## RECURSIVE CALL TO CURRENT SCRIPT
			## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			## WARNING: start the tgen_mutation_burden script [we call the current script actually so we MUST NOT provide again the option --pqarm otherwise infinite loop]
			mycmd="bash $(readlink -en $0) -t ${THREADS} --vcf ${VCF} --nbam ${NBAM} --tbam ${TBAM} --outprefix ${PREFIX_PQARM} --region ${COORD} --region-name ${PQARM_COORD_NAME} ${LOPT[@]}"
			echo "${mycmd}"
			eval ${mycmd}
			ev=$?
			if [[ ${ev} -eq 0 ]] ;
			then
				LOSF="${LOSF[@]} ${OSF}"
			else
				echo -e "ERROR: tgen_mutation_burden script FAILED for coordinates ${COORD} with bam pair << ${NBAM}:::${TBAM} >>; " 1>&2 ;
				echo ${OSF} >> ${LOG_PQARM_FAILURE}
			fi
		done < ${PQARMS_FILE}
	else
		echo -e "using GNU parallel ..."
		type parallel >/dev/null 2>&1 || { echo >&2 "Require  \"Parallel\"   but it's not in the PATH. Please install GNU parallel or add it to your PATH before restarting the jobs.  Aborting."; exit 1; }
		if [[ ${PADDED_BED} != "" ]] ; then LOPT="--bed ${PADDED_BED}" ; fi
		if [[ ${FORCE} == "yes" ]] ; then LOPT="${LOPT[@]} --force " ; fi
		# echo -e "unbuffer parallel --xapply -j ${GNU_PARALLEL}  bash $(readlink -en $0) -t ${THREADS} --library \"${LIBRARY}\" --sample \"${SAMPLE}\" --rg \"${READ_GROUP}\" --vcf ${VCF} --nbam ${NBAM} --tbam ${TBAM} --outprefix parallel__${PREFIX}pqarms.{2} --region {1} --region-name {3} ${LOPT[@]} ::: `cat ${PQARMS_FILE} | cut -f1 | tr '\n' ' ' ` ::: `cat ${PQARMS_FILE} | cut -f1 | sed 's/[:-]/_/g ' | tr '\n' ' '` ::: `cat ${PQARMS_FILE} | cut -f2 | sed 's/ \+//g; s/[:-]/_/g' | tr '\n' ' '` "
		if [[ ${READ_GROUP} == "" ]] ; then READ_GROUP="NA" ; fi
		if [[ ${LIBRARY} == "" ]] ; then LIBRARY="NA" ; fi
		if [[ ${SAMPLE} == "" ]] ; then SAMPLE="NA" ; fi
		unbuffer parallel --xapply -j ${GNU_PARALLEL}  bash $(readlink -en $0) -t ${THREADS}  --vcf ${VCF} --nbam ${NBAM} --tbam ${TBAM} --outprefix parallel__${PREFIX}pqarms.{2} --region {1} --region-name {3} --verbose ${LOPT[@]} --library "${LIBRARY}" --sample "${SAMPLE}" --rg "${READ_GROUP}" ::: `cat ${PQARMS_FILE} | cut -f1 ` ::: `cat ${PQARMS_FILE} | cut -f1 | sed 's/[:-]/_/g ' ` ::: `cat ${PQARMS_FILE} | cut -f2 | sed 's/ \+//g; s/[:-]/_/g'`  &>log_Parallel_MutBurden.${PREFIX}${TBAM_BN}.txt

		if [[ $? -eq 0 ]] ;
		then
			echo -e "creating LOSF ..."
			LOSF=()
			LOSF=( $(find * -mindepth 0 -maxdepth 0 -type f -name "parallel__${PREFIX}pqarms*${TBAM_BN}.coverage.mutation_burden.stats.txt" | tr '\n' ' ' | sort -k1V ) )
			if [[ $? -ne 0 || ${#LOSF[@]} -eq 0 ]] ; then echo -e "\nERROR: No FILE have been Created; Check your inputs; Aborting." ; fi
			echo -e "LOSF ===> \n ${LOSF[@]} \n <=== LOSF"
		else
			echo -e "\nERROR: section << GNU parallel >> in tgen_mutation_burden script FAILED with bam pair << ${NBAM}:::${TBAM_BN} >>; \nAborting without having created the ggplot\n " | tee -a  >> log_Parallel_MutBurden.${PREFIX}${TBAM_BN}.txt
			exit 1
		fi
	fi

	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## loop for processing each region or interval is DONE, now we do the plot
	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if [[ "${VERBOSE}" == "yes" ]] ; then echo "List_File(s) Run with GNU_PARALLEL:  ${LOSF[@]}" ; echo -e "FAILURE for Region File(s) if exists == '${LFAILURE[@]}' " ; fi
	if [[ ${SKIP_PLOT} == "no" ]]
	then
		echo -e "\nMaking PLOT ...\n"
		## concatenating the parallel data and creating the file for plotting mutation burden using R script
		CONCAT_PQARMS_OUTFILE="${PREFIX}${TBAM_BN}.coverage.mutation_burden.summary.stats.for_Rplot.txt"
		if [[ ${VERBOSE} == "yes" ]];
		then
			echo -e "VERBOSE: CONCAT_PQARMS_OUTFILE == ${CONCAT_PQARMS_OUTFILE}" ;
			echo -e "cat ${LOSF[@]} | grep -vE \"^#\" | grep -A1 MutationsCount | grep -vE \"MutationsCount|\-\-\" | sort -k10V | awk '{ split(\$10,C,\":\") ; \$1=C[1]\"\t\"\$10\"\t\"\$1 ; OFS=\"\t\" ; print }' | cut -f1-5,13-14 | sed 's/\t/\t/' > ${CONCAT_PQARMS_OUTFILE}"
		fi

		cat ${LOSF[@]} | grep -vE "^#" | grep -A1 MutationsCount | grep -vE "MutationsCount|\-\-" | sort -k10V | awk '{ split($10,C,":") ; $1=C[1]"\t"$10"\t"$1 ; OFS="\t" ; print }' | cut -f1-5,13-14 | sed 's/\t/\t/' > ${CONCAT_PQARMS_OUTFILE}
		echo -e "chrom\tlocus\tMutationsCount\tCoverage\tMutationBurdenPerMillionBase\tregion_name\tpct_dbsnp" > header_for__coverage.mutation_burden.summary.stats.txt_file__.txt
		## running R script to plot mutation burden per arm
		LOPTforRscript=" --infile ${CONCAT_PQARMS_OUTFILE} "  ## LOPT stands for listOPTion
		if [[ ${ADD_DBSNP_TO_PLOT} == "yes" ]] ; then LOPTforRscript="${LOPTforRscript[@]} --addPctDbsnp " ; fi
		if [[ $( head -n1 ${CONCAT_PQARMS_OUTFILE} | cut -f1 | grep -cE "^chr")  -eq 0 ]] ; then LOPTforRscript="${LOPTforRscript[@]} --noPrefixChr" ; fi

		## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		## running script to make plot for each region given
		## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		echo -e "\nRscript $(dirname $(readlink -en $0))/mutation_burden_per_arm.R ${LOPTforRscript[@]}\n"
		Rscript $(dirname $(readlink -en $0))/mutation_burden_per_arm.R ${LOPTforRscript[@]}
		if [[ $? -ne 0 ]] ; then echo -e "\nERROR: ggplot creation FAILED ; Check log file to get the Rscript command in order to rerun it manually, after fixing why it failed; Aborting." ; exit 1 ; fi
	fi

	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## FINALIZING OUTPUT from PARALLEL RUNS
	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	SQUISH_OUTFILE=${PREFIX}${TBAM_BN}.coverage.mutation_burden.stats.txt
	if  [[ ${GNU_PARALLEL} -lt 2 ]]
	then
			if [[ ${OUTFILE} != "" ]] ;
			then
				echo -e "renaming output file ... "
				echo -e "mv \${SQUISH_OUTFILE} \${OUTFILE} <==>  mv ${SQUISH_OUTFILE} ${OUTFILE}"
				mv ${SQUISH_OUTFILE} ${OUTFILE}
				echo -e "writing results into file ... \t${OUTFILE}"
				pipelinize "${OUTFILE}" "${PIPELINE}"
			fi
			echo -e "writing results into file ... \t${SQUISH_OUTFILE}"
			pipelinize "${SQUISH_OUTFILE}" "${PIPELINE}"
	else
		if [[ ${OUTFILE} == "" ]] ; then OUTFILE=${SQUISH_OUTFILE} ; fi
		echo -e "writing results into file ... \t${OUTFILE}"
		true >"${OUTFILE}"
		## add headers to the output statistics summary file aka OUTFILE_STATS;
		if [[ ${PIPELINE} == "yes" ]]
		then
			echo -e "# ${ORIGINAL_COMMAND}" >> "${OUTFILE}"
			echo -e "##MAXIMUM_GIVEN_CALLABLE_POSITIONS=${MAX_CALLABLE_SPACE}" >> "${OUTFILE}"
			echo -e "##FILTER=mindepth[FS+RS]>=${MINDEPTH} && minReadPerStrandInSTumorOnly[TFS+TRS]>=${MINREADSPERSTRAND} ; TFS=TumorForwardStrand;TRS=TumorReverseStrand" >> "${OUTFILE}"
			echo -e "# Started on:  $( date )" >> "${OUTFILE}"
		else
			echo -e "##Started on:  $( date )" >> "${OUTFILE}"
			echo -e "##${ORIGINAL_COMMAND}" >> "${OUTFILE}"
			echo -e "##MAXIMUM_GIVEN_CALLABLE_POSITIONS=${MAX_CALLABLE_SPACE}" >> "${OUTFILE}"
			echo -e "##FILTER=mindepth[FS+RS]>=${MINDEPTH} && minReadPerStrandInSTumorOnly[TFS+TRS]>=${MINREADSPERSTRAND} ; TFS=TumorForwardStrand;TRS=TumorReverseStrand" >> "${OUTFILE}"
		fi
#		head -n7 $(find -maxdepth 1 -type f -name "parallel__*${TBAM}.coverage.mutation_burden.stats.txt" | sort | head -n 1) ;
		(  echo -e "\n## METRICS CLASS\t" ; (echo -e "$(cat parallel__*${TBAM}.coverage.mutation_burden.stats.txt  | grep -vE "^#|READ_GROUP" | sed '/^$/d' | sort -k10V | awk '{FS=OFS="\t" ; SUM1+=$1 ; SUM2+=$2 ; SUM3+=$3 ; C4[$4]=$4 } END { FS=OFS="\t";  printf SUM1"\t"SUM2"\t"SUM1/SUM2*1000000} ')\t${MINDEPTH}\t${MINREADSPERSTRAND}\t${MIN_BASE_QUAL}\t${MIN_MAP_QUAL}\t${NBAM_BN}\t${TBAM_BN}\tNA\tNA\tNA\t${SAMPLE}\t${LIBRARY}\t${READ_GROUP}"  ) ) >> ${OUTFILE}
		pipelinize "${OUTFILE}" "${PIPELINE}"
	fi

	exit
fi
## @@@@@@@@@@@@@@@@@@@@@@@@@
## END PARALLEL RUNS
## @@@@@@@@@@@@@@@@@@@@@@@@@

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## SERIAL RUN SECTION
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## 4 cases:  1) No Region  and No BED given ->ERROR ; 2) No REGION but BED given ; 3) REGION given AND BED given ; 4) REGION given but NO BED given: this case is used either for debugging or automatically when '--pqarms' option used; The parallel cmd calls that case;
## When '--pqarms' file is given, the regions are the p and q arms of each chromosome; This is A SPECIAL CASE we process above using multiple cpus (gnu-parallel) or not (gnu-paralle -lt 2)
## Setting variables with options to be used with "samtools view" and "samtools depth" depends upon whether REGION is provided or not
if [[ ${REGION_INTERVAL} == "" ]] ;
then
	if [[ "${PADDED_BED}" == "" ]] ;
	then
		echo -e "\nERROR: FNF ; BED File MUST be provided using option --bed , particularly if no region has been provided using --region; Aborting." ;
		exit 1  ;
	fi
fi


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## WHAT REGIONS ?
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if [[ ${REGION_INTERVAL} != "" && "${PADDED_BED}" != "" ]] ;
then
	echo -e "## Both a Region and a Bed file have been provided" ;
	echo -e "## We use the Region because '--region' has precedence over '--bed'";
	echo -e "## If you want to use the regions in the BED file only, do not use the '--region' option"
	echo -e "--> Processing Region now ..."
	echo -e "## Using the provided region ${REGION_INTERVAL} ; Therefore Samtools depth should run faster.\nsamtools depth in progress ..."
	CONTIG=$(echo ${REGION_INTERVAL} | cut -d":" -f1)
	START=$(echo ${REGION_INTERVAL} | cut -d":" -f2 | cut -d"-" -f1)
	STOP=$(echo ${REGION_INTERVAL} | cut -d":" -f2 | cut -d"-" -f2)
	MAX_CALLABLE_SPACE=$( echo "${STOP}-${START}" | bc -l | awk '{printf "%0.0f", $0}' )
	samtools_view_options="-h -q ${MIN_MAP_QUAL} -M " ;
	samtools_depth_options="-q ${MIN_BASE_QUAL} " ;
	PADDED_BED=bed_region_$$_${TBAM_BN}_temp.bed
	echo -e "${CONTIG}\t${START}\t${STOP}" > "${PADDED_BED}"
elif [[ ${REGION_INTERVAL} != "" && "${PADDED_BED}" == "" ]]
then
	echo -e "Only Region has been provided with '--region' option" ;
	echo -e "## Using the provided region ${REGION_INTERVAL} ;  No BED file provided ; Therefore Samtools depth should run faster.\nsamtools depth in progress ..."
	echo -e "## or can be b/c we use parallel and made a recursive call to the script. See Parallel Section of the script $0"
	CONTIG=$(echo ${REGION_INTERVAL} | cut -d":" -f1)
	START=$(echo ${REGION_INTERVAL} | cut -d":" -f2 | cut -d"-" -f1)
	STOP=$(echo ${REGION_INTERVAL} | cut -d":" -f2 | cut -d"-" -f2)
	MAX_CALLABLE_SPACE=$( echo "${STOP}-${START}" | bc -l | awk '{printf "%0.0f", $0}' )
	samtools_view_options="-h -q ${MIN_MAP_QUAL} -M" ;
	samtools_depth_options="-q ${MIN_BASE_QUAL} " ;

	echo -e "we convert the region or given interval into a temp bed file"
	PADDED_BED=bed_region_$$_${TBAM_BN}_temp.bed
	echo -e "${CONTIG}\t${START}\t${STOP}" > "${PADDED_BED}"
elif [[ "${PADDED_BED}" != "" ]]
then
	echo -e "## only BED file provided with '--bed'"
	echo -e "NOTE: you provided a bed file: make sure the contig names match between the bam files and the bed file."
	MAX_CALLABLE_SPACE=$( cat ${PADDED_BED} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' )
	echo -e "samtools depth in progress ..."
	samtools_view_options="-h -q ${MIN_MAP_QUAL} -L ${PADDED_BED} -M"
	samtools_depth_options="-b ${PADDED_BED} -q ${MIN_BASE_QUAL}"
else
  echo -e "something went wrong. Should not have happened; Contact your Admin or the Tool's Author" ;
  exit 1
fi

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MISC PROCESSING STEPS and CHECKS
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "MAX_CALLABLE_SPACE=${MAX_CALLABLE_SPACE}"
if [[ ${MAX_CALLABLE_SPACE} -eq 0 ]] ; then echo -e "ERROR: \$MAX_CALLABLE_SPACE equals ZERO ; Something is wrong; Check inputs; Aborting now"; exit 1 ; fi

## we want to print the position with zero coverage as well in the unfiltered compressed log file so we need to use the '-aa' option in samtools depth
if [[ ${MINDEPTH} -eq 0 ]] ; then samtools_depth_options="${samtools_depth_options} -aa " ; fi

## if provided, adding THREADS to samtools view options
if [[ ${THREADS} -gt 1 ]] ; then samtools_view_options=""${samtools_view_options[@]}" -@ ${THREADS}" ; fi

## For DEBUGGING or logging
if [[ "${VERBOSE}" == "yes" ]] ;
then
	echo -e "\nsamtools_depth_options == \"${samtools_depth_options}\""
	echo -e "samtools_view_options == \"${samtools_view_options}\""

	echo -e "samtools depth ${samtools_depth_options} " ;
	echo -e "samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL}" ;
	echo -e "samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL}" ;
	echo -e "samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL}" ;
	echo -e "samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL}" ;
	echo -e "awk -F\"\\\t\" -v MINDEPTH=${MINDEPTH} -v MINREADSPERSTRAND=${MINREADSPERSTRAND} '{ OFS=\"\\\t\" ; NORMAL_SUM=\$3+\$4 ; TFS=\$5 ; TRS=\$6 ; TUMOR_SUM=TFS+TRS ; if(NORMAL_SUM>=MINDEPTH && TUMOR_SUM>=MINDEPTH && TFS>=MINREADSPERSTRAND && TRS>=MINREADSPERSTRAND ) { print \$0 } }' | wc -l \n"
fi

## Test here if user specified an outfilename
if [[ ${OUTFILE} != "" ]] ; then OUTFILE_STATS=${OUTFILE} ; fi


## (re-)init stats outfile if any previous version existed before running samtools depth
## This is a feature that can be used separately from the TMB calculation
## (even though the TMB calculation is pretty fast) --> Therefore I would advice using the
## option '--keep-unflt-outfile' instead.
if [[ ${RUN_ONLY_DEPTH_CAPTURE} == "no"  ]]
then
	true >"${OUTFILE_STATS}"
	## getting the read coverage over the PADDED_BED file which represents the exome targets + padded target + expanded up to exons
	echo -e "## Samtools depth in progress ..."
	## add headers to the output statistics summary file aka OUTFILE_STATS;
	if [[ ${PIPELINE} == "yes" ]]
	then
		echo -e "# $(echo -e "$(realpath $0)") $@" >> "${OUTFILE_STATS}"
		echo -e "##MAXIMUM_GIVEN_CALLABLE_POSITIONS=${MAX_CALLABLE_SPACE}" >> "${OUTFILE_STATS}"
		echo -e "##FILTER=minMapQual=${MIN_MAP_QUAL} && minBaseQual=${MIN_BASE_QUAL}" >> "${OUTFILE_STATS}"
		echo -e "##FILTER=mindepth[FS+RS]>=${MINDEPTH} && minReadPerStrandInSTumorOnly[TFS+TRS]>=${MINREADSPERSTRAND} ; TFS=TumorForwardStrand;TRS=TumorReverseStrand" >> "${OUTFILE_STATS}"
		echo -e "# Started on:  $( date )" >> "${OUTFILE_STATS}"
  else
		echo -e "##$(echo -e "$(realpath $0)") $@" >> "${OUTFILE_STATS}"
		echo -e "##MAXIMUM_GIVEN_CALLABLE_POSITIONS=${MAX_CALLABLE_SPACE}" >> "${OUTFILE_STATS}"
		echo -e "##FILTER=minMapQual=${MIN_MAP_QUAL} && minBaseQual=${MIN_BASE_QUAL}" >> "${OUTFILE_STATS}"
		echo -e "##FILTER=mindepth[FS+RS]>=${MINDEPTH} && minReadPerStrandInSTumorOnly[TFS+TRS]>=${MINREADSPERSTRAND} ; TFS=TumorForwardStrand;TRS=TumorReverseStrand" >> "${OUTFILE_STATS}"
		echo -e "##Started on:  $( date )" >> "${OUTFILE_STATS}"
	fi
fi





## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MAIN SAMTOOLS COMMAND
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if [[ ${KEEP_OUTFILE_UNFILTCOUNT_FSRS} != "yes" ]] ;
then
	coverage=$( ( samtools depth ${samtools_depth_options} \
	<(samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL} 2>err.F0x10_F0x400.Constitutional.$$.log ) \
	<(samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL} 2>err.f0x10_F0x400.Constitutional.$$.log ) \
	<(samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL} 2>err.F0x10_F0x400.Case.$$.log ) \
	<(samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL} 2>err.f0x10_F0x400.Case.$$.log ) ) \
	| awk -F"\t" -v MINDEPTH=${MINDEPTH} -v MINREADSPERSTRAND=${MINREADSPERSTRAND} '{ OFS="\t" ; NORMAL_SUM=$3+$4 ; TFS=$5 ; TRS=$6 ; TUMOR_SUM=TFS+TRS ; if(NORMAL_SUM>=MINDEPTH && TUMOR_SUM>=MINDEPTH && TFS>=MINREADSPERSTRAND && TRS>=MINREADSPERSTRAND ) { print $0 } }' | wc -l )
	ev=$?
else
	echo -e "NOTE: you asked to keep the intermediate UNfiltered depth values."
	echo -e "NOTE: columns order in << *.CallablePositions.unflt.txt.gz >> file:"
	echo -e "NOTE: order is: CONTIG\tPOSITION\tNormalForwardStrand\tNormalReverseStrand\tTumorForwardStrand\tTumorReverseStrand"
	coverage=$( ( samtools depth ${samtools_depth_options} \
	<(samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL} 2>err.F0x10_F0x400.Constitutional.$$.log ) \
	<(samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${NBAM} ${REGION_INTERVAL} 2>err.f0x10_F0x400.Constitutional.$$.log ) \
	<(samtools view -F0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL} 2>err.F0x10_F0x400.Case.$$.log ) \
	<(samtools view -f0x10 -F0x400 ${samtools_view_options[@]} ${TBAM} ${REGION_INTERVAL} 2>err.f0x10_F0x400.Case.$$.log ) ) \
	| tee >(gzip -c - > ${OUTFILE_UNFILTCOUNT_FSRS} ) | awk -F"\t" -v MINDEPTH=${MINDEPTH} -v MINREADSPERSTRAND=${MINREADSPERSTRAND} '{ OFS="\t" ; NORMAL_SUM=$3+$4 ; TFS=$5 ; TRS=$6 ; TUMOR_SUM=TFS+TRS ; if(NORMAL_SUM>=MINDEPTH && TUMOR_SUM>=MINDEPTH && TFS>=MINREADSPERSTRAND && TRS>=MINREADSPERSTRAND ) { print $0 } }' | wc -l )
	ev=$?
fi
echo -e "NOTE: About Filtering the positions to calculate the Mutation Burden; we only keep the positions that are true for the following statement:"
echo -e "NOTE: (NORMAL_SUM>=MINDEPTH && TUMOR_SUM>=MINDEPTH && TFS>=MINREADSPERSTRAND && TRS>=MINREADSPERSTRAND )"
for F in ${LIST_ERROR_FILES[@]}  ;
do
	if [[ $(cat ${F} | wc -l ) -ne 0 ]] ;
		then echo -e "*** WARNING--WARNING *** ---> a warning or error with << samtools view >> got raised; check file << ${F} >> ; \nNOTE: it might just be because the index bai file is older than the data bam itself; but if not, you MUST determine if warning or error might have impacted the mutation burden calculation by giving incorrect results" ;
		else rm ${F} ;
	fi ;
done

if [[ ${ev} -ne 0 ]] ;
then
	echo -e "\nERROR: samtools depth FAILED; ev = ${ev} ;  Aborting! ; " ; echo -e "FAILED with code ${ev}" >> "${OUTFILE_STATS}" ; exit 1 ;
fi
echo -e "## samtools depth ... DONE.\n"

if [[ ${RUN_ONLY_DEPTH_CAPTURE} == "yes" ]] ;
then
	echo -e "You asked to run only the coverage capture using 'samtools depth' command(s) and Exit."  ;
	echo -e "Depth for your Regions are available here: ${PWD} "
	exit ;
fi




## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## CALCULATIONS PER INTERVAL
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "## Calculation Section ... "
## grabbing number of mutations in VCF
## first Checking if vcf is a block compressed vcf?
MADE_GZ="no"
if [[ ${VCF##*.} != "gz" ]] ;
then
	echo -e "bzipping and tabixing vcf ..."
	bgzip -c ${VCF} > ${VCF}_tmp.$$.vcf.gz ;
	tabix -p vcf ${VCF}_tmp.$$.vcf.gz
	VCF=${VCF}_tmp.$$.vcf.gz
	MADE_GZ="yes"
fi


echo -e "## Step RSID COUNT ... with FLAG_DBSNP = ${FLAG_DBSNP}"
if [[ ${REGION_INTERVAL} == ""  ]] ;
then
	echo -e "No other regions of interest were given beside the (padded)_BED file; So we proceed with all the Regions given in BED file ${PADDED_BED}"
	region="all_regions_in_target_bed"
	##On the Debate side ...
	##Important note worth mentioned:
	##this cmd will count the lines in vcf; if the same position exist on TWO lines in VCF because Two different mutations have been found at that position, it will be counted as two different mutations (use sort -u in the command below to count it as one)
	##if the Two mutations are mentioned on hte same line but separated with a comma in the ALT column, these two mutations will be counted as one, since we count the number of lines (modify the command to split and count the number of ALT alleleas as the number of Mutations)

	mutations=$( bedtools intersect -a <( zcat -f ${VCF} | grep -vE "^#" | awk '{FS=OFS="\t" ; $2=$2-1"\t"$2 ; print $0 }' ) -b ${PADDED_BED} | cut -f1-3,5-6 | wc -l )

	if [[ "${FLAG_DBSNP}" == "RS" ]]
	then
		echo -e "\${FLAG_DBSNP} == \"RS\" "
		RSID_COUNT=$(echo "$(zcat -f ${VCF} | grep -vE "^#" | cut -f3 | grep -E "^rs[0-9]*" | wc -l )" )
	else
		echo -e "\${FLAG_DBSNP} == \"${FLAG_DBSNP}\" "
		RSID_COUNT=$(echo "$(zcat -f ${VCF} | grep -vE "^#" | cut -f8 | grep -E "${FLAG_DBSNP}" | wc -l )" )
	fi

else
	echo -e "(padded)_BED and/or Region(s) were given ... "
	echo -e "We capture the number of mutations within the Region(s) of Interest ... "
	region=${REGION_INTERVAL}

	REGION_BEDFILE=temp_region_interval_$$_$(echo -e "${REGION_INTERVAL}" | sed 's/[:-]/_/g').bed
	echo -e "${REGION_INTERVAL}" | sed 's/[:-]/\t/g' > ${REGION_BEDFILE}

	##See above the 'On the Debate Side' commented out paragraph as well for additional information about how the mutations are counted
	mutations=$( bedtools intersect -a <( zcat -f ${VCF} | grep -vE "^#" | awk '{FS=OFS="\t" ; $2=$2-1"\t"$2 ; print $0 }' ) -b ${PADDED_BED} | bedtools intersect -a stdin -b ${REGION_BEDFILE} | cut -f1-3,5-6 | wc -l )


	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	##  CAPTURE COUNT of DATA for adding it to output stats file
	## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if [[ $? -eq 0 ]] ; then rm ${REGION_BEDFILE} &>/dev/null ; fi
	if [[ "${FLAG_DBSNP}" == "RS" ]]
	then
		echo -e "value for \${FLAG_DBSNP} == \"RS\" "
		RSID_COUNT=$(echo "$(tabix ${VCF} ${REGION_INTERVAL} | grep -vE "^#" | cut -f3 | grep -E "^rs[0-9]*" | wc -l )" )
	else
		echo -e "value for \${FLAG_DBSNP} == \"${FLAG_DBSNP}\" "
		RSID_COUNT=$(echo "$(tabix ${VCF} ${REGION_INTERVAL} | grep -vE "^#" | cut -f8 | grep -E "${FLAG_DBSNP}" | wc -l )" )
	fi
fi  # end else of if REGION_INTERVAL


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## CALCULATE % DBSNP for STATS OUTFILE if '--pqarms' option got used
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "RSID_COUNT == ${RSID_COUNT}"
echo -e "Checking whether mutations number is zero or more for calculating the percentage of DBSNP in that region ... "
if [[ ${mutations} -gt 0 && ${RSID_COUNT} -gt 0 ]] ;
then
	echo -e "number of mutations found is: ${mutations}"
	PCT_DBSNP=$( echo "${RSID_COUNT}/${mutations}*100" | bc -ql | awk '{printf "%f", $0}' )
else
	echo -e "number of mutations found is Zero or RSID is Zero [ nmut: ${mutations} :: rsidCount: ${RSID_COUNT} ]"
	PCT_DBSNP="NA"
fi

echo -e "PCT_DBSNP == ${PCT_DBSNP}"

if [[ ${MADE_GZ} == "yes"  && -e ${VCF}_tmp.$$.vcf.gz ]] ;
then
	echo -e "removing temp gz file and their index..."
	rm -r ${VCF}_tmp.$$.vcf.gz ${VCF}_tmp.$$.vcf.gz.tbi ;
	echo -e "ev rm: $? ; " ;
fi


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## MUTATION BURDEN CALCULATION
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "calculate Mutation Burden per million bases ..."
echo -e "coverage: ${coverage}"
echo -e "mutations: ${mutations}"

if [[ ${coverage} != "" && ${coverage} -gt 0 && ${mutations} -gt 0 ]] ;
then
	mutBurden=$(echo "(${mutations}/${coverage})*1000000" | bc -l | awk '{printf "%f", $0}' )
else
	echo -e "\nWARNING: Either Number of mutations is zero -- or --  Unknown Coverage value -- or -- Coverage value is Zero --> Division par Zero! ; Check stats file for more details; NOTE: that might be a true results actually\n" 1>&2
	mutBurden=NA
fi


if [[ ${READ_GROUP} == "NA" ]] ; then READ_GROUP="" ; fi
if [[ ${LIBRARY} == "NA" ]] ; then LIBRARY="" ; fi
if [[ ${SAMPLE} == "NA" ]] ; then SAMPLE="" ; fi

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## output statistics to OUTFILE_STATS file the same formatting way Picard
## does to ease the future parsing of the file with the existing
## picard's wrapper used in TGen's pipeline (JK's request)
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
echo -e "writing results into file ... \t${OUTFILE_STATS}"
echo -e "\n## METRICS CLASS\t
MutationsCount\tCoverage\tMutationBurdenPerMillionBase\tMinReadDepth\tTumorMinReadDepthPerStrand\tMinBaseQual\tMinMapQual\tNBAM\tTBAM\tregion\tregion_name\tpct_dbsnp\tSAMPLE\tLIBRARY\tREAD_GROUP
${mutations}\t${coverage}\t${mutBurden}\t${MINDEPTH}\t${MINREADSPERSTRAND}\t${MIN_BASE_QUAL}\t${MIN_MAP_QUAL}\t${NBAM_BN}\t${TBAM_BN}\t${region}\t${REGION_NAME}\t${PCT_DBSNP}\t${SAMPLE}\t${LIBRARY}\t${READ_GROUP}" >> ${OUTFILE_STATS}

pipelinize "${OUTFILE_STATS}" ${PIPELINE}

## print summary if verbose enabled ; for developper, useful when debugging
if [[ "${VERBOSE}" == "yes" ]] ; then echo -e "\n##Summary print ...\n" ; cat ${OUTFILE_STATS} ; fi
echo -e "Results in:\t${PWD}"
echo $(date)
