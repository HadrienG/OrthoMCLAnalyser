#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Written by Hadrien Gourlé. Feel free to use and modify"""

import os
import argparse
import numpy as np
import itertools
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import GC
import sys # testing


def parser():
	Description = "Orthomcl output parser"
	parser = argparse.ArgumentParser(description=Description)
	parser.add_argument(
		"--input",
		"-i",
		help="Input file",
		metavar="[groups.txt]",
		required=True
	)
	parser.add_argument(
		"--gb_dir",
		"-g",
		help="Directory containing the .gb of the genomes \
		used to build the cluster",
		metavar="[GbDir]",
		required=True
	)
	parser.add_argument(
		"--genome_list",
		"-l",
		help="GenomeList.txt",
		metavar="[file]",
		required=True
	)
	parser.add_argument(
		"--outdir",
		"-o",
		default="OutDir",
		help="Output Directory (default: %(default)s)"
	)
	args = parser.parse_args()
	return args


class Genome(object):
	"""class representing a genome coming from a .gb file"""
	def __init__(self, element):
		super(Genome, self).__init__()
		self.type = ""
		self.filename = element[0]
		self.name = element[0].split("/")[-1].split(".")[0]
		self.code = element[1]
		self.seq = ""
		self.size = 0
		self.coding_size = 0
		self.CDSList = []
		self.TranslatedCDS = {}  # TranslatedCDS[locus_tag] = DNAsequence
		if self.filename.endswith((".gb",".gbk",".gbf")):
			self.type = "genbank"
			self.file = SeqIO.parse(open(self.filename), "genbank")
			print "Parsing %s" % (self.filename)
			for record in self.file:
				self.seq += str(record.seq)
				self.size += (len(record.seq))
				for feature in record.features:
					if feature.type == "CDS":
						# We add 1 to the aa length to take into account the stop codon
						self.coding_size += (len("".join(feature.qualifiers["translation"])) + 1) * 3
						self.CDSList.append(feature)
						self.TranslatedCDS["".join(feature.qualifiers["locus_tag"])] = feature.location.extract(record).seq
			self.GC = GC(self.seq)

		elif self.filename.endswith((".fasta",".faa")):
			self.type = "fasta"
			self.file = SeqIO.parse(open(self.filename), "fasta")
			print "Parsing %s" % (self.filename)
			accession = None
			for record in self.file:
				header, description, sequence = record.id, record.description, str(record.seq)
				if accession is None:
					accession = header.split("|")[1].split(".")[0]
				self.coding_size += (len("".join(sequence)) + 1) * 3

			for record in SeqIO.parse(Entrez.efetch("nucleotide", id=accession, rettype="gbwithparts"),"genbank"):
				self.seq += str(record.seq)
				self.size += (len(record.seq))
				for feature in record.features:
					if feature.type == "CDS":
						self.CDSList.append(feature)
						self.TranslatedCDS["".join(feature.qualifiers["locus_tag"])] = feature.location.extract(record).seq
			self.GC = GC(self.seq)


class Protein(object):
	"""class representing a protein, member of a Cluster."""
	# Protein.strain		= The 3 or 4 char code used by Orthomcl
	def __init__(self, header, args, Genomes):
		super(Protein, self).__init__()
		self._info = header
		self.strain = header.split("|")[0]
		self._genbank = self._load_file(args, header, Genomes)
		self.feature = self._genbank[2]
		self.product = "".join(self._genbank[0])
		self.sequence = "".join(self._genbank[1])
		self.length = len(self.sequence)

	def _load_file(self, args, header, Genomes):
		Genome = (x for x in Genomes if x.code == header.split("|")[0])
		for genome in Genome:
			if genome.type == "fasta":
				for feature in genome.CDSList:
					if feature.qualifiers["protein_id"][0].split(".")[0] == header.split("|")[1].split(".")[1][7:]:
						if "product" in feature.qualifiers:
							product = feature.qualifiers["product"]
							sequence = feature.qualifiers["translation"]
							feat = feature
						else:
							product = ["unknown"]
							sequence = feature.qualifiers["translation"]
							feat = feature
						break
			elif genome.type == "genbank":
				for feature in genome.CDSList:
					if feature.qualifiers["locus_tag"][0] == header.split("|")[1]:
						if "product" in feature.qualifiers:
							product = feature.qualifiers["product"]
							sequence = feature.qualifiers["translation"]
							feat = feature
						else:
							product = ["unknown"]
							sequence = feature.qualifiers["translation"]
							feat = feature
						break
		return product, sequence, feat


class Cluster(object):
	"""class representing a cluster of proteins"""
	def __init__(self, line, args, Genomes):
		self.id = line.split()[0][:-1]
		self.protlist = line.split()[1:]
		self.proteins = []
		self.strainlist = set()

		for prot in self.protlist:
			x = Protein(prot, args, Genomes)
			print "\tProtein with header %s: done." % (x._info)
			self.proteins.append(x)

		for protein in self.proteins:
			self.strainlist.add(protein.strain)

		if len(self.strainlist) < len(Genomes):
			self.type = "PAN"
		else:
			self.type = "CORE"

		# Some useful stats about the cluster
		self.lenlist = [x.length for x in self.proteins]

		self.maxlength = max(self.lenlist)
		self.minlength = min(self.lenlist)
		self.mean = np.mean(self.lenlist)
		self.stdev = np.std(self.lenlist)

		# Representative protein of the cluster: the longest or the first one
		self.repr = "".join(next(
			x.sequence for x in self.proteins if x.length == self.maxlength)
		)


def extract_clusters(args, Genomes):
	# return a list of Cluster objects from the groups.txt file
	Clusters = []
	with open(args.input) as groups:
		for line in groups:
			print "Extracting info for Cluster"
			x = Cluster(line, args, Genomes)
			print "Extracting info for Cluster %s: done." % (x.id)
			Clusters.append(x)
	return Clusters


def identify_singletons(Clusters, genome):
	# Some proteins are not present in the orthoMCL clusters. This function
	# aims to identify them. They should be included in each strain pan-genome.
	CompleteSet = set(genome.CDSList)
	PartialProteinList = [x for x in list(
		itertools.chain(*[x.proteins for x in Clusters])) if x.strain == genome.code]
	PartialFeatureSet = set(protein.feature for protein in PartialProteinList)

	Singletons = CompleteSet - PartialFeatureSet
	return Singletons  # Return a set of features


def cluster2fasta(cluster, output):
	# Writes a .faa containing all the members of a cluster
	Output = open(output, "w")
	for name, prot in zip(cluster.protlist, cluster.proteins):
		Output.write(">%s\n%s\n" % (name, prot.sequence))
	Output.close()
	print "%s.faa written" % (cluster.id)


def cluster_summary(cluster, output):
	# Writes a file summarizing the content of a cluster
	form = """ID: %s\nTYPE: %s\nMEMBERS: %s\nNUMBER OF MEMBERS: %s\n
	Mean Length: %.2f (std dev: %.2f)
	Min Length: %s
	Max Length: %s\n\nREPRESENTANT OF THE CLUSTER:\n%s\n\nPRODUCTS:\n%s
	""" % (
		cluster.id,
		cluster.type,
		" ".join(cluster.strainlist),
		len(cluster.proteins),
		cluster.mean,
		cluster.stdev,
		cluster.minlength,
		cluster.maxlength,
		cluster.repr,
		"".join(["\t" + x.product + "\n" for x in cluster.proteins])
	)
	Output = open(output, "w")
	Output.write(form)
	Output.close()
	print "Summary for cluster %s written" % (cluster.id)


def overall_summary(args, Clusters):
	form = """Summary of the Clusters content:\n
	Genomes used: %s
	Number of Clusters: %i
	"Core" Clusters: %i
	"Pan" Clusters: %i\n
	Total number of proteins in Clusters: %i
	Number of "Core" proteins: %i
	Number of "Pan" proteins: %i\n
	List of "Core" Clusters: %s
	List of "Pan" Clusters: %s
	""" % (
		" ".join(os.listdir(args.gb_dir)),
		len(Clusters),
		len([x.id for x in Clusters if x.type == "CORE"]),
		len([x.id for x in Clusters if x.type == "PAN"]),
		len(list(itertools.chain(*[x.proteins for x in Clusters]))),
		len(list(itertools.chain(*[x.proteins for x in Clusters
											if x.type == "CORE"]))),
		len(list(itertools.chain(*[x.proteins for x in Clusters
											if x.type == "PAN"]))),
		[x.id for x in Clusters if x.type == "CORE"],
		[x.id for x in Clusters if x.type == "PAN"]
	)
	Output = open(args.outdir + "Clusterlist.txt", "w")
	Output.write(form)
	Output.close()
	print "%sClusterlist.txt written." % (args.outdir),


def get_core_genome(args, Clusters, genome):
	CoreFeatures = [x.feature for x in list(
		itertools.chain(*[x.proteins for x in Clusters if x.type == "CORE"]))
		if x.strain == genome.code]  # List of features

	CoreProteins = [feature.qualifiers["translation"] for feature in CoreFeatures]
	CoreDNA = {"".join(feature.qualifiers["locus_tag"]): genome.TranslatedCDS["".join(feature.qualifiers["locus_tag"])]
		for feature in CoreFeatures}
	MergedCoreDNA = "".join(str(seq) for seq in CoreDNA.values())
	CoreGC = GC(MergedCoreDNA)

	# Writing .fna and .faa for core genome in the Genomes dir
	OutCoreProt = open(args.outdir + "/Genomes/" + genome.name + ".Core.faa", "w")
	OutCoreDNA = open(args.outdir + "/Genomes/" + genome.name + ".Core.fna", "w")
	OutMergedDNA = open(args.outdir + "/Genomes/" + genome.name + ".CoreMerged.fna", "w")

	for feature in CoreFeatures:
		OutCoreProt.write(">%s\n%s\n" % ("".join(feature.qualifiers["locus_tag"]), feature.qualifiers["translation"][0]))
	OutCoreProt.close()
	for key, value in CoreDNA.iteritems():
		OutCoreDNA.write(">%s\n%s\n" % (key, value))
	OutCoreDNA.close()
	OutMergedDNA.write(">%s\n%s\n" % (genome.name + "_CORE", MergedCoreDNA))
	OutMergedDNA.close()

	# TODO WRITE A README

	return {"CoreDNA":CoreDNA, "MergedCoreDNA":MergedCoreDNA, "CoreGC":CoreGC, "CoreProteins":CoreProteins}


def get_pan_genome(args, Clusters, genome):
	PanFeatures = [x.feature for x in list(
		itertools.chain(*[x.proteins for x in Clusters if x.type == "PAN"]))
		if x.strain == genome.code] + [feature for feature in identify_singletons(Clusters, genome)]
		# List of features

	PanProteins = [feature.qualifiers["translation"] for feature in PanFeatures]
	PanDNA = {"".join(feature.qualifiers["locus_tag"]): genome.TranslatedCDS["".join(feature.qualifiers["locus_tag"])]
		for feature in PanFeatures}
	MergedPanDNA = "".join(str(seq) for seq in PanDNA.values())
	PanGC = GC(MergedPanDNA)

	# Writing .fna and .faa for core genome in the Genomes dir
	OutPanProt = open(args.outdir + "/Genomes/" + genome.name + ".Pan.faa", "w")
	OutPanDNA = open(args.outdir + "/Genomes/" + genome.name + ".Pan.fna", "w")
	OutMergedDNA = open(args.outdir + "/Genomes/" + genome.name + ".PanMerged.fna", "w")

	for feature in PanFeatures:
		OutPanProt.write(">%s\n%s\n" % ("".join(feature.qualifiers["locus_tag"]), feature.qualifiers["translation"][0]))
	OutPanProt.close()
	for key, value in PanDNA.iteritems():
		OutPanDNA.write(">%s\n%s\n" % (key, value))
	OutPanDNA.close()
	OutMergedDNA.write(">%s\n%s\n" % (genome.name + "_PAN", MergedPanDNA))
	OutMergedDNA.close()

	return {"PanDNA":PanDNA, "MergedPanDNA":MergedPanDNA, "PanGC":PanGC, "PanProteins":PanProteins}


def per_genome_summary(args, Clusters, Genomes):
	form = "Genome\tN. of proteins\tN. of Clusters\tN. of Singletons\n"
	form2 = "Genome\tGC\tSize\tCodingSize\t%Coding\tCoreSize\tPanSize\t%Core\t%Pan\tGC Core\tGC Pan\n"
	Output = open(args.outdir + "Summary.txt", "w")
	for genome in Genomes:
		form += "%s\t%i\t%i\t%s\n" % (
			genome.code,
			len([x for x in list(itertools.chain(*[x.proteins for x in Clusters]))
										if x.strain == genome.code]),
			len([x.id for x in Clusters if genome.code in x.strainlist]),
			len(identify_singletons(Clusters, genome))
		)

		CoreGenomeDic = get_core_genome(args, Clusters, genome)
		PanGenomeDic = get_pan_genome(args, Clusters, genome)
		# {"PanDNA":PanDNA, "MergedPanDNA":MergedPanDNA, "PanGC":PanGC, "PanProteins":PanProteins}

		form2 += "%s\t%.2f\t%i\t%i\t%.2f\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\n" % (
								genome.code,
								genome.GC,
								genome.size,
								genome.coding_size,
								genome.coding_size / float(genome.size) * 100,
								len(CoreGenomeDic["MergedCoreDNA"]),
								len(PanGenomeDic["MergedPanDNA"]),
								len(CoreGenomeDic["MergedCoreDNA"]) / float(genome.coding_size) * 100,
								len(PanGenomeDic["MergedPanDNA"]) / float(genome.coding_size) * 100,
								CoreGenomeDic["CoreGC"],
								PanGenomeDic["PanGC"]
		)
	Output.write(form2)
	Output.close()
	print "%sSummary.txt written.\n" % (args.outdir)
	print form, "\n", form2, "\n"




def makedir(args):
	# Creates output directory
	try:
		os.makedirs(args.outdir + "/CoreGenome")
		os.makedirs(args.outdir + "/PanGenome")
		os.makedirs(args.outdir + "/Genomes")
		print "Output Directory Created"
	except OSError:
		if not os.path.isdir(args.outdir):
			raise


def main():
	args = parser()
	makedir(args)
	Entrez.email ="mail@example.com"

	GenomeList=[]
	with open(args.genome_list,"r") as file:
		for line in file:
			GenomeList.append((args.gb_dir+line.split(",")[0], line.split(",")[2]))
	# List of objects Genome
	Genomes = [Genome(element) for element in GenomeList]
	# print Genomes

	# List of objects Cluster
	Clusters = extract_clusters(args, Genomes)

	# Writing of the output
	for cluster in Clusters:
		if cluster.type == "PAN":
			cluster_summary(cluster, args.outdir + "/PanGenome/" + cluster.id + ".txt")
			cluster2fasta(cluster, args.outdir + "/PanGenome/" + cluster.id + ".faa")
		else:
			cluster_summary(cluster, args.outdir + "/CoreGenome/" + cluster.id + ".txt")
			cluster2fasta(cluster, args.outdir + "/CoreGenome/" + cluster.id + ".faa")
	overall_summary(args, Clusters)
	print

	# Summary of the analysis + Core/Pan size calculations
	per_genome_summary(args, Clusters, Genomes)


if __name__ == '__main__':
	main()
