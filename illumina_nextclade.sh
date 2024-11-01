
OUTDIR=${1}

mkdir ${OUTDIR}/results/Nextclade/

grep -v sylvatic ${OUTDIR}/results/virus_calls.tsv | sed '1d' | cut -f 1,2,4 > ${OUTDIR}/.serotypes.tsv

while read SAMPLE CONSENSUS SEROTYPE
do

    nextclade run \
		--input-dataset "$(dirname "$(realpath "$0")")/DENV_Nextclade/${SEROTYPE}" \
		--output-all ${OUTDIR}/results/Nextclade/${SAMPLE} \
		${OUTDIR}/results/consensus_sequences/${CONSENSUS}

done < ${OUTDIR}/.serotypes.tsv

rm -f ${OUTDIR}/.serotypes.tsv

(cat $(ls ${OUTDIR}/results/Nextclade/*/nextclade.tsv | head -n 1) && tail -n +2 -q ${OUTDIR}/results/Nextclade/*/nextclade.tsv) > ${OUTDIR}/results/Nextclade/nextclade.tsv