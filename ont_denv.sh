
################################################################################################
###### NOTE: This script is intended to be run within an activated ARTIC conda environment #####
################################################################################################
########## MODIFIED FROM https://github.com/josephfauver/DENV_MinION_Script by CALUM ###########
############################ IF IT DOESN'T WORK, BLAME HIM #####################################
################################################################################################
######################### SPECIFY VARIABLES NEEDED TO RUN THE PIPELINE #########################
################################################################################################

## assign the appropriate medaka model based on sequencing and basecalling parameters
MEDAKA_MODEL=r1041_e82_400bps_sup_v4.2.0
# capute working directory
BASE_DIR="$(pwd)"

################################################################################################
######################### SHOULD NOT NEED TO EDIT CODE BELOW THIS LINE #########################
################################################################################################

## positional argument specifying absolute path to directory in which fastqs are located:
## this should be a directory containing one file per sample, (e.g. "DMG2404857_2.fastq.gz")
READS_DIR=${BASE_DIR}/$1
## positional argument specifying absolute path to output directory:
OUTDIR=${BASE_DIR}/$2
## path to scheme directory:
SCHEME_DIR="$(dirname "$(realpath "$0")")/DENV_Schemes"

# make output directory
mkdir ${OUTDIR}
# make fofn of input FASTQs
ls ${READS_DIR}/*.fastq.gz | sed 's,.*/,, ; s,.fastq.gz,,' > ${OUTDIR}/.temp_fofn
######################################################
#### reports how many samples are being processed ####
######################################################
echo ""
NUM_SAMPLES=$(cat ${OUTDIR}/temp_fofn | wc -l)
echo "***********************************************"
echo "** TOTAL NUMBER OF SAMPLES TO PROCESS: " ${NUM_SAMPLES}
echo "***********************************************"
echo ""
#################################################################
#### runs the pipeline, looping over each line in input list ####
#################################################################
cd ${OUTDIR}

while read SAMPLE
do

	echo "***********************************************"
	echo "** PROCESSING" ${SAMPLE}
	echo "***********************************************"
	echo ""
	echo "***********************************************"
	echo "** DETERMINING SEROTYPE FOR" ${SAMPLE}"..."
	echo "***********************************************"
	echo ""

	DENV1_REF=${SCHEME_DIR}/DENV1/V1/DENV1.reference.fasta
	DENV2_REF=${SCHEME_DIR}/DENV2/V1/DENV2.reference.fasta
	DENV3_REF=${SCHEME_DIR}/DENV3/V1/DENV3.reference.fasta
	DENV4_REF=${SCHEME_DIR}/DENV4/V1/DENV4.reference.fasta

	minimap2 -ax asm20 ${DENV1_REF} ${READS_DIR}/${SAMPLE}.fastq.gz | samtools sort -o ${SAMPLE}_DENV1_sorted.bam -T ${SAMPLE}_DENV1.tmp
	minimap2 -ax asm20 ${DENV2_REF} ${READS_DIR}/${SAMPLE}.fastq.gz | samtools sort -o ${SAMPLE}_DENV2_sorted.bam -T ${SAMPLE}_DENV2.tmp
	minimap2 -ax asm20 ${DENV3_REF} ${READS_DIR}/${SAMPLE}.fastq.gz | samtools sort -o ${SAMPLE}_DENV3_sorted.bam -T ${SAMPLE}_DENV3.tmp
	minimap2 -ax asm20 ${DENV4_REF} ${READS_DIR}/${SAMPLE}.fastq.gz | samtools sort -o ${SAMPLE}_DENV4_sorted.bam -T ${SAMPLE}_DENV4.tmp
	
	DENV1=$(samtools view -c -F 260 ${SAMPLE}_DENV1_sorted.bam)
	echo ""
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV1:" ${DENV1}
	echo "***********************************************"
	echo ""

	DENV2=$(samtools view -c -F 260 ${SAMPLE}_DENV2_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV2:" ${DENV2}
	echo "***********************************************"
	echo ""

	DENV3=$(samtools view -c -F 260 ${SAMPLE}_DENV3_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV3:" ${DENV3}
	echo "***********************************************"
	echo ""

	DENV4=$(samtools view -c -F 260 ${SAMPLE}_DENV4_sorted.bam)
	echo "***********************************************"
	echo "** NO. READS FOR" ${SAMPLE}" MAPPED TO DENV4:" ${DENV4}
	echo "***********************************************"
	echo ""
	
	MAX=$DENV1
	MAX_VAR="DENV1"
	if (( $DENV2 > $MAX )) ; then MAX=$DENV2 ; MAX_VAR="DENV2" ; fi
	if (( $DENV3 > $MAX )) ; then MAX=$DENV3 ; MAX_VAR="DENV3" ; fi
	if (( $DENV4 > $MAX )) ; then MAX=$DENV4 ; MAX_VAR="DENV4" ; fi
	SEROTYPE=${MAX_VAR}

	if (( ${MAX} == 0 ))
	then
		echo "********************************************************"
		echo "** NO READS MAPPED; ASSINGING SEROTYPE DENV1 BY DEFAULT."
		echo "********************************************************"
		echo ""
		SEROTYPE="DENV1"
	else
		
		echo "      ***********************************"
		echo "      **" ${SAMPLE} "IS SEROTYPE" ${SEROTYPE}
		echo "      ***********************************"
		echo ""
	fi

	rm *.bam

	mkdir ${SAMPLE}
	cd ${SAMPLE}

	artic minion \
		--scheme-directory ${SCHEME_DIR} \
		--read-file ${READS_DIR}/${SAMPLE}.fastq.gz \
		--normalise 400 \
		--medaka \
		--medaka-model ${MEDAKA_MODEL} \
		${SEROTYPE}/V1 ${SAMPLE}

	echo ""
	echo "***********************************************************************"
	echo "** ANALYSIS FOR SAMPLE" ${SAMPLE} "FINISHED. CHECK ARTIC OUTPUT FOR ERRORS."
	echo "***********************************************************************"
	echo ""
	echo ""

	nextclade run \
		--input-dataset "$(dirname "$(realpath "$0")")/DENV_Nextclade/${SEROTYPE}" \
		--output-all Nextclade \
		${SAMPLE}.consensus.fasta

	cat ${SAMPLE}.coverage_mask.txt.1.depths ${SAMPLE}.coverage_mask.txt.2.depths > ${SAMPLE}.depths.txt
	
	csvtk plot line \
		-t \
		-H \
		-x 3 \
		-y 4 \
		-g 2 \
		--point-size 0.01 \
		--legend-top \
		--xlab Position \
		--ylab Depth \
		-o ${SAMPLE}.depths.pdf \
		${SAMPLE}.depths.txt

	cd ..

done < .temp_fofn

(cat $(ls */Nextclade/*.tsv | head -n 1) && tail -n +2 -q */Nextclade/*.tsv) > nextclade.tsv

rm -f temp_fofn 

echo "ALL ANALYSES COMPLETE. GOODBYE."
echo""