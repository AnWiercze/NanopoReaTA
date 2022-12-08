function getMappedReads {
	bam="$1/*bam"
	outFile=$2
	echo "Sample\tnum_reads\tnum_mapped_reads" > $outFile
	for i in $bam; do
		numReads=$( samtools view -c $bam )
		numMappedReads=$( samtools view -c -F 260 $bam )
		sample=$( basename $bam ) 
		sample=${sample/.bam/}
		echo "$sample\t$numReads\t$numMappedReads" >> $outFile 
	done
}
echo $1 
echo $2
getMappedReads $1 $2/mapping_stats.txt
