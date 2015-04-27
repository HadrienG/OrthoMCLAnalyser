#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO


def parser():
    Description = "Genbank to .faa converter"
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument(
                        "--input",
                        "-i",
                        help="input file",
                        metavar="[file.gb]",
                        required=True
    )
    parser.add_argument(
                        "--output",
                        "-o",
                        help="output prefix",
                        metavar="[file]",
                        required=True
    )
    args = parser.parse_args()
    return args


def GenbankToFaa(args):
    Input = open(args.input, "r")
    Output = open("%s.faa" % args.output, "w")

    for Seq in SeqIO.parse(Input, "genbank"):
        print "Dealing with GenBank record %s" % Seq.id
        for SeqFeature in Seq.features:
            if SeqFeature.type == "CDS":
                assert len(SeqFeature.qualifiers['translation']) == 1
                if "locus_tag" in SeqFeature.qualifiers:
                    Output.write(">%s|%s\n%s\n" % (
                                                        SeqFeature.qualifiers
                                                        ['locus_tag'][0],
                                                        Seq.name,
                                                        SeqFeature.qualifiers
                                                        ['translation'][0]))

    Input.close()
    Output.close()


def main():
    args = parser()
    GenbankToFaa(args)


if __name__ == "__main__":
    main()
