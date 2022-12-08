#### All Functions from RSeQC package version 4.0.0
# Some functions got minimal changes
import sys
import json
import argparse
import math

# Argument
parser = argparse.ArgumentParser()
parser.add_argument('--bed')
parser.add_argument('--output_dir')
args = parser.parse_args()

def percentile_list(N):
    """
    Find the percentile of a list of values.
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @return - the list of percentile of the values
    """
    if not N:return None
    if len(N) <100: return N
    per_list=[]
    for i in range(1,101):
        k = (len(N)-1) * i/100.0
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            per_list.append( int(N[int(k)])  )
        else:
            d0 = N[int(f)] * (c-k)
            d1 = N[int(c)] * (k-f)
            per_list.append(int(round(d0+d1)))
    return per_list


def genebody_percentile(refbed, outDir, mRNA_len_cut = 100):
	
	g_percentiles = {}
	transcript_count = 0
	for line in open(refbed,'r'):
		try:
			if line.startswith(('#','track','browser')):continue  
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			geneID = geneName
				
			exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
			exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
			exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
			exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
			transcript_count += 1
		except:
			print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
			continue
		gene_all_base=[]
		mRNA_len =0
		flag=0
		for st,end in zip(exon_starts,exon_ends):
			gene_all_base.extend(list(range(st+1,end+1)))		#1-based coordinates on genome
		if len(gene_all_base) < mRNA_len_cut:
			continue
		g_percentiles[geneID] = (chrom, strand, percentile_list(gene_all_base))	#get 100 points from each gene's coordinates
	
	with open(outDir + '/g_percentiles.json', 'w') as fp:
		json.dump(g_percentiles, fp)

## Bed file from https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE.v38.bed.gz + gunzip
#refbed = "/path/to/bed/file/hg38_GENCODE.v38.bed"
genebody_percentile(args.bed, args.output_dir)
