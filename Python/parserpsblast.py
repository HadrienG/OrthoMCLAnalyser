#!/usr/bin/env python

import glob
from Bio.Blast import NCBIXML


E_VALUE_THRESH = 0.00001

CoreList = glob.glob("*Core*")
PanList = glob.glob("*Pan*")

Output = open("../CoreCogList.txt","w")
for xml in CoreList:
    with open(xml,"r") as file:
        for record in NCBIXML.parse(file):
            if record.alignments:
                # print "QUERY: %s..." % record.query[:60]
                for align in record.alignments :
                    Output.write(align.hit_def.split(",")[0]+"\n")
                    break
                    # for hsp in align.hsps :
                    #     print " %s HSP, e=%f, from position %i to %i" \
                    #     % (align.hit_def.split(",")[0], hsp.expect, hsp.query_start, hsp.query_end)
                    #     assert hsp.expect <= E_VALUE_THRESH
Output.close()

Output = open("../PanCogList.txt","w")
for xml in PanList:
    with open(xml,"r") as file:
        for record in NCBIXML.parse(file):
            if record.alignments:
                for align in record.alignments :
                    Output.write(align.hit_def.split(",")[0]+"\n")
                    break
Output.close()
