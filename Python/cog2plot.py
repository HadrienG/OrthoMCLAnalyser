#!/usr/bin/env python

import re


def ConstructFunDic(Fun, CogNames):
    Codes = {fun.split("\t")[0]: fun.strip().split("\t")[1] for fun in Fun}
    Dic = {}
    CogNames.readline()
    for cog in CogNames:
        line = cog.split("\t")
        if len(line[1]) == 1:
            Dic[line[0]] = (line[1], Codes[line[1]])
        elif re.match("^[A-Za-z]*$", line[1]):
            for code in line[1]:
                Dic[line[0]] = (code, Codes[code])
    return Dic

def main():
    Fun = open("CogInfo/fun2003-2014.tab","r")
    CogNames = open("CogInfo/cognames2003-2014.tab","r")
    FunDic = ConstructFunDic(Fun, CogNames)

    with open("PanCoglist.txt", "r") as PanCog, open("CoreCoglist.txt", "r") as CoreCog:
        OutPan = open("PanCog.tmp.txt", "w")
        OutCore = open("CoreCog.tmp.txt", "w")

        for line in CoreCog:
            outline = "%s\t%s\t%s\n" % (
                line.strip(),
                FunDic[line.strip()][0],
                FunDic[line.strip()][1]
                )
            OutCore.write(outline)

        for line in PanCog:
            outline = "%s\t%s\t%s\n" % (
                line.strip(),
                FunDic[line.strip()][0],
                FunDic[line.strip()][1]
                )
            OutPan.write(outline)

        OutPan.close()
        OutCore.close()

if __name__ == "__main__":
    main()
