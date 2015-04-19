#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Written by Hadrien Gourlé. Feel free to use and modify"""

import os
import argparse
import numpy as np
import itertools
from Bio import SeqIO


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
		"--genome_number",
		"-n",
		help="Number of genomes used to build the clusters",
		type=int,
		metavar="[int]",
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
	def __init__(self, genbank):
		super(Genome, self).__init__()
		self.file = genbank
		self.name = genbank.split("/")[-1].split(".")[0]
		self.size = 0
		self.coding_size = 0
		self.CDSList = []

		Genbank = SeqIO.parse(open(self.file), "genbank")
		for record in Genbank:
			self.size += (len(record.seq))
			for feature in record.features:
				if feature.type == "CDS":
					self.coding_size += len("".join(feature.qualifiers["translation"])) * 3
					self.CDSList.append(feature)


class Protein(object):
	"""class representing a protein, member of a Cluster."""
	# Contains:
	# Protein.strain		= The 3/4 chars code used by Orthomcl
	# Protein.product		= The function of the protein
	# Protein.sequence		= The sequence of the protein
	def __init__(self, header, args, Genomes):
		super(Protein, self).__init__()
		self._info = header
		self.strain = header.split("|")[0]
		self._genbank = self._load_genbank(args, header, Genomes)
		self.feature = self._genbank[2]
		self.product = "".join(self._genbank[0])
		self.sequence = "".join(self._genbank[1])
		self.length = len(self.sequence)

	def _load_genbank(self, args, header, Genomes):
		Genome = (x for x in Genomes if x.name == header.split("|")[0])
		for genome in Genome:
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

		if len(self.strainlist) < args.genome_number:
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
		itertools.chain(*[x.proteins for x in Clusters])) if x.strain == genome.name]
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
	form = """Summary of the analysis:\n
	Genomes used: %s
	Number of Clusters: %i
	"Core" Clusters: %i
	"Pan" Clusters: %i\n
	Total number of proteins: %i
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

	print form


def per_genome_summary(args, Clusters, Genomes):
	# Output = open(args.outdir + "/clusterlist.txt", "w")
	form = "Genome\tN. of proteins\tN. of Clusters\tN. of Singletons\n"
	# GenomeList = os.listdir(args.gb_dir)
	for genome in Genomes:
		form += "%s\t%i\t%i\t%s\n" % (
			genome.name,
			len([x for x in list(itertools.chain(*[x.proteins for x in Clusters]))
										if x.strain == genome.name]),
			len([x.id for x in Clusters if genome.name in x.strainlist]),
			len(identify_singletons(Clusters, genome))
			# [x.id for x in Clusters if genome.strip(".gb") in x.strainlist]
		)
	# Output.write(form)
	# Output.close()
	print form


def makedir(args):
	# Creates output directory
	try:
		os.makedirs(args.outdir + "/CoreGenome")
		os.makedirs(args.outdir + "/PanGenome")
		print "Output Directory Created"
	except OSError:
		if not os.path.isdir(args.outdir):
			raise


def main():
	args = parser()
	makedir(args)

	# List of objects Genome
	Genomes = [Genome(args.gb_dir + gb) for gb in os.listdir(args.gb_dir)]

	# List of objects Cluster
	Clusters = extract_clusters(args, Genomes)

	for cluster in Clusters:
		if cluster.type == "PAN":
			cluster_summary(cluster, args.outdir + "/PanGenome/" + cluster.id + ".txt")
			cluster2fasta(cluster, args.outdir + "/PanGenome/" + cluster.id + ".faa")
		else:
			cluster_summary(cluster, args.outdir + "/CoreGenome/" + cluster.id + ".txt")
			cluster2fasta(cluster, args.outdir + "/CoreGenome/" + cluster.id + ".faa")

	print
	per_genome_summary(args, Clusters, Genomes)

	# overall_summary(args, Clusters)
	print
	# size_calculations(Clusters)

	# GENOMES
	print "Genome\tSize\tCodingSize\t%Coding\tCoreSize\tPanSize"
	# CoreSize + PanSize < Coding size. TODO retrieve the singletons
	for genome in Genomes:
		print "%s:\t%i\t%i\t%f\t%i\t%i" % (
								genome.name,
								genome.size,
								genome.coding_size,
								genome.coding_size / float(genome.size) * 100,
								sum([x.length for x in list(
									itertools.chain(*[x.proteins for x in Clusters if x.type == "CORE"]))
									if x.strain == genome.name]) * 3,
								sum([x.length for x in list(
									itertools.chain(*[x.proteins for x in Clusters if x.type == "PAN"]))
									if x.strain == genome.name]) * 3
		)

	# singleton identification DONE.
	# TODO calculate the length of the singletons and add it to PANSize


if __name__ == '__main__':
	main()
