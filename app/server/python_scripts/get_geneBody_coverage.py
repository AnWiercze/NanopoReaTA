import sys
from numpy import std,mean
import operator
import pysam
import collections
import pandas as pd
import json
from os.path import basename
import glob

import argparse

parser = argparse.ArgumentParser(description='Run Gene Body Coverage.')
parser.add_argument('--bamList', nargs='+')
parser.add_argument('--gene')
parser.add_argument('--converted_gtf')
parser.add_argument('--output_dir')

args = parser.parse_args()

print(args)

def genebody_coverage(bam, position_list):
	'''
	position_list is dict returned from genebody_percentile
	position is 1-based genome coordinate
	'''
	samfile = pysam.Samfile(bam, "rb")
	chr_bool = False
	for read in samfile.fetch():
		if "chr" in str(read):
			chr_bool = True
			break
	aggreagated_cvg = collections.defaultdict(int)
	
	gene_finished = 0
	for chrom, strand, positions in list(position_list.values()):
		chrom = chrom.replace("chr","")
		if chr_bool == True:
			chrom = f"chr{chrom}"
		coverage = {}
		for i in positions:
			coverage[i] = 0.0
		chrom_start = positions[0]-1
		if chrom_start <0: chrom_start=0
		chrom_end = positions[-1]
		print(chrom_start,chrom_end)
		
		for pileupcolumn in samfile.pileup(chrom, chrom_start, chrom_end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in positions:
				continue
			if pileupcolumn.n == 0:
				coverage[ref_pos] = 0
				continue				
			cover_read = 0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				if pileupread.alignment.is_duplicate:continue
				cover_read +=1
			coverage[ref_pos] = cover_read
		tmp = [coverage[k] for k in sorted(coverage)]
		if strand == '-':
			tmp = tmp[::-1]
		for i in range(0,len(tmp)):
			aggreagated_cvg[i] += tmp[i]
		gene_finished += 1
		
		if gene_finished % 100 == 0:
			print("\t%d transcripts finished\r" % (gene_finished), end=' ', file=sys.stderr)
	return 	aggreagated_cvg

def reduceDict(gp, transOfInterest):
	gp_sub = dict((k, gp[k]) for k in transOfInterest['transcript_id'].values
           if k in gp)
	return	gp_sub

def getGeneCoverage(bamFilesList, geneOfInterest, gtfFile, outDir):
	# Output file name
	outFile = outDir + "/" + "samples.geneBodyCoverage.txt"
	OUT1 = open(outFile	,'w')
	print("Percentile\t" + '\t'.join([str(i) for i in range(1,101)]), file=OUT1)
	# Load gtf to map transcripts to genes
	gtf = pd.read_csv(gtfFile)
	# Load dict containing transcript percentiles 
	with open(outDir + "/g_percentiles.json") as json_file:
		gp = json.load(json_file)
	#for i in gp.keys():
	#	gp[i][0] = gp[i][0].replace("chr", "")
	#print(gtf.head())
	gtf_sub = gtf[['gene_id', 'transcript_id', 'transcript_name', 'gene_name']].drop_duplicates()
	transOfInterest = gtf_sub[gtf_sub['gene_id']== geneOfInterest].dropna()
	gp_sub = reduceDict(gp, transOfInterest)
	file_container = []
	for bamfile in bamFilesList:
		print(bamfile)
		cvg = genebody_coverage(bamfile, gp_sub)
		print(cvg)
		if len(cvg) == 0:
			print("\nCannot get coverage signal from " + basename(bamfile) + ' ! Skip', file=sys.stderr)
			continue
		tmp = valid_name(basename(bamfile).replace('.bam',''))	# scrutinize R identifer
		if file_container.count(tmp) == 0:
			print(tmp + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		else:
			print(tmp + '.' + str(file_container.count(tmp)) + '\t' + '\t'.join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
		file_container.append(tmp)
	OUT1.close()

def valid_name(s):
	'''make sure the string 's' is valid name for R variable'''
	symbols = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.'
	digit = '0123456789'
	rid = '_'.join(i for i in s.split())	#replace space(s) with '_'
	if rid[0] in digit:rid = 'V' + rid
	tmp = ''
	for i in rid:
		if i in symbols:
			tmp = tmp + i
		else:
			tmp = tmp + '_'
	return tmp

bamFilesList = glob.glob(str(args.bamList[0]) + "*.bam")
print(bamFilesList)
geneOfInterest = args.gene
gtf = args.converted_gtf
outDir = args.output_dir
#bamFilesList = ['/media/anna/MinION_Drive/run_dir_22_09/bam_genome_merged/ERR6053055.bam', '/media/anna/MinION_Drive/run_dir_22_09/bam_genome_merged/ERR6053056.bam', '/media/anna/MinION_Drive/run_dir_22_09/bam_genome_merged/ERR6053097.bam', '/media/anna/MinION_Drive/run_dir_22_09/bam_genome_merged/ERR6053098.bam']
#geneOfInterest = 'ENSG00000160182.3'
#gtf='/media/anna/MinION_Drive/run_dir_22_09/converted_gtf.csv'
#outDir='/media/anna/MinION_Drive/run_dir_22_09/'

getGeneCoverage(bamFilesList, geneOfInterest, gtf, outDir)
	
