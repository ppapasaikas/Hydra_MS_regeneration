$NTHREAD=16;
$INDEXDIR="Hydra_vulgaris_r102/starIndex_Hydra_vulgaris_ncbi_r102/";
$GTF="genomes/Hydra_vulgaris_r102/Hydra_vulgaris_ncbi_r102_processed.gtf";
$trimCommand="--clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
$tpm="Basic";
$readCommand="--readFilesCommand zcat";
$FASTQDIR="my_FASTQ_directory";
$OUTDIR="my_OUTDIR";
    

opendir(DIR, $FASTQDIR) or die $!;
my @FQfiles = grep {
/^*.gz/			 	

 && -f "$FASTQDIR/$_"   # and is a file
} readdir(DIR);


foreach $fq (@FQfiles){
next if $fq=~/Undeterm/;
$c++;

$prefix=$1 if  $fq=~/(.+).fastq.gz/;
$fqpath=$FASTQDIR . $fq;
$starout=$OUTDIR . $prefix;

#Skip if bam output exists:
next if -f $OUTDIR . 'BAM/' . $prefix . 'Aligned.sortedByCoord.out.bam';

print "STAR --runThreadN $NTHREAD --gensomeDir $INDEXDIR --sjdbGTFfile $GTF --readFilesIn $fqpath $readCommand $trimCommand --outSJfilterReads Unique --outFilterType BySJout --outFilterMultimapNmax 5 --alignSJoverhangMin 8 --alignSJDBoverhangMin 4 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --seedSearchStartLmax 50 --twopassMode $tpm --outFileNamePrefix $starout --quantMode GeneCounts\n";
print `STAR --runThreadN $NTHREAD --genomeDir $INDEXDIR --sjdbGTFfile $GTF --readFilesIn $fqpath $readCommand $trimCommand --outSJfilterReads Unique --outFilterType BySJout --outFilterMultimapNmax 5 --alignSJoverhangMin 8 --alignSJDBoverhangMin 4 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --seedSearchStartLmax 50 --twopassMode $tpm --outFileNamePrefix $starout --quantMode GeneCounts --sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFtagExonParentGene gene_id --genomeChrBinNbits 12 --genomeSAsparseD 2`;
}