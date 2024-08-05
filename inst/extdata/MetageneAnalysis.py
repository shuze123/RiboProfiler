import collections

##############################################################################################################
# for mapping to genome qc analysis
##############################################################################################################

def MetageneAnalysis(geneInfo, inputFile, outputFile, mode="st", type="codon", cdslength=600, expression=75, exclude=90):
    ################################################################
    # 1. choose which mode to analysis
    ################################################################
    if mode == "st":
        mylist = range(-50, 1501)
        region = [i for i in range(-51, 1501, 3)]
        codonPos = [i for i in range(-17, 500)]
    elif mode == "sp":
        mylist = range(-1500, 51)
        region = [i for i in range(-1501, 51, 3)]
        codonPos = [i for i in range(-500, 17)]
    else:
        print("please give the st/sp mode")
        return

    # save info
    rangeDict = collections.defaultdict(float)
    countDict = collections.defaultdict(float)

    ################################################################
    # 2. extract gene feature info
    ################################################################
    featureLenDict = {}

    # open gene info file
    with open(geneInfo, "r") as input:
        for line in input:
            _, _, gene_id, trans_id, _, _, _, _, utr5, cds, utr3 = line.split()
            utr5, cds, utr3 = int(utr5), int(cds), int(utr3)
            cdsStart, cdsEnd = (utr5 + 1), (utr5 + cds - 2)
            featureLenDict[trans_id] = (cdsStart, cdsEnd, cds)

    ################################################################
    # 3. filter gene CDS length and expression higher than threshold
    ################################################################
    gene_infoDict = collections.defaultdict(float)
    filtedGeneDict = {}

    # open gene info file
    with open(inputFile, "r") as input:
        for line in input:
            fileds = line.split()
            trans_id = fileds[2]
            transpos = int(fileds[3])
            rpm = float(fileds[4])

            # filter CDS > 600 nt gene
            if trans_id in featureLenDict:
                cdsStart, cdsEnd, cdsLength = featureLenDict[trans_id]
                if cdsLength >= cdslength and cdsLength % 3 == 0:
                    key = f"{trans_id}:{cdsLength - exclude}"
                    if mode == "st":
                        if cdsStart + exclude <= transpos <= cdsEnd:
                            gene_infoDict[key] += rpm
                    else:
                        if cdsStart <= transpos <= cdsEnd - exclude:
                            gene_infoDict[key] += rpm

    # filter CDS expression >= 75
    for key, val in gene_infoDict.items():
        if val >= expression:
            meanNorm = val / int(key.split(":")[1])
            filtedGeneDict[key] = meanNorm

    ################################################################
    # 4. Meta-gene analysis from start codon
    ################################################################
    with open(inputFile, "r") as input:
        for line in input:
            fileds = line.split()
            trans_id = fileds[2]
            transpos = int(fileds[3])
            density = float(fileds[4])
            cdsStart, cdsEnd, cdsLength = featureLenDict[trans_id]
            id = f"{trans_id}:{cdsLength - exclude}"

            if id in filtedGeneDict:
                if mode == "st":
                    reldist = transpos - cdsStart
                elif mode == "sp":
                    reldist = transpos - cdsEnd

                if mylist[0] <= reldist <= mylist[-1]:
                    reads = density / filtedGeneDict[id]
                    rangeDict[reldist] += reads
                    countDict[reldist] += 1

    ################################################################
    # 5. output data
    ################################################################
    fullDict = {}
    for k in rangeDict:
        fullDict[k] = rangeDict[k] / countDict[k] if countDict[k] != 0 else 0

    meanDensity = sum(fullDict.values()) / len(fullDict)
    new_fullDict = {k: v / meanDensity for k, v in fullDict.items()}

    if type == "codon":
        codonDict = {}
        for i in range(len(region) - 1):
            codonRegion = range(region[i], region[i] + 3)
            count = 0
            codonDensity = 0
            for j in codonRegion:
                if j in new_fullDict:
                    codonDensity += new_fullDict[j]
                    count += 1
            codonDict[codonPos[i]] = codonDensity / count if count != 0 else 0
        finalDict = codonDict
    else:
        finalDict = new_fullDict

    with open(outputFile, "w") as outFileP:
        for pos, meanDensity in sorted(finalDict.items()):
            outFileP.write(f"{pos}\t{meanDensity}\n")


##############################################################################################################
# for mapping to transcriptome qc analysis
##############################################################################################################

def MetageneAnalysis_ontrans(inputFile, outputFile, mode="st", type="codon", cdslength=600, expression=75, exclude=90):
    if mode == "st":
        mylist = range(-50, 1501)
        region = [i for i in range(-51, 1501, 3)]
        codonPos = [i for i in range(-17, 500)]
    elif mode == "sp":
        mylist = range(-1500, 51)
        region = [i for i in range(-1501, 51, 3)]
        codonPos = [i for i in range(-500, 17)]
    else:
        print("please give the st/sp mode")
        return

    rangeDict = collections.defaultdict(float)
    countDict = collections.defaultdict(float)

    gene_infoDict = collections.defaultdict(float)
    filtedGeneDict = {}

    with open(inputFile, "r") as input:
        for line in input:
            fileds = line.split()
            gene_name, _, trans_id, cdsStart, cdsEnd, _ = fileds[1].split("|")
            pos = int(fileds[2])
            cdsLength = int(cdsEnd) - int(cdsStart)
            density = float(fileds[4])

            if cdsLength > cdslength and cdsLength % 3 == 0:
                key = f"{trans_id}:{cdsLength - exclude}"
                if mode == "st":
                    if int(cdsStart) + exclude <= pos <= int(cdsEnd):
                        gene_infoDict[key] += density
                else:
                    if int(cdsStart) <= pos <= int(cdsEnd) - exclude:
                        gene_infoDict[key] += density

    for key, val in gene_infoDict.items():
        if val > expression:
            meanNorm = val / int(key.split(":")[1])
            filtedGeneDict[key] = meanNorm

    with open(inputFile, "r") as input:
        for line in input:
            fileds = line.split()
            gene_name, _, trans_id, cdsStart, cdsEnd, _ = fileds[1].split("|")
            pos = int(fileds[2])
            cdsLength = int(cdsEnd) - int(cdsStart)
            density = float(fileds[4])
            id = f"{trans_id}:{cdsLength - 90}"

            if id in filtedGeneDict:
                if mode == "st":
                    reldist = pos - int(cdsStart)
                elif mode == "sp":
                    reldist = pos - int(cdsEnd)

                if mylist[0] <= reldist <= mylist[-1]:
                    reads = density / filtedGeneDict[id]
                    rangeDict[reldist] += reads
                    countDict[reldist] += 1

    fullDict = {}
    for k in rangeDict:
        fullDict[k] = rangeDict[k] / countDict[k] if countDict[k] != 0 else 0

    meanDensity = sum(fullDict.values()) / len(fullDict)
    new_fullDict = {k: v / meanDensity for k, v in fullDict.items()}

    if type == "codon":
        codonDict = {}
        for i in range(len(region) - 1):
            codonRegion = range(region[i], region[i] + 3)
            count = 0
            codonDensity = 0
            for j in codonRegion:
                if j in new_fullDict:
                    codonDensity += new_fullDict[j]
                    count += 1
            codonDict[codonPos[i]] = codonDensity / count if count != 0 else 0
        finalDict = codonDict
    else:
        finalDict = new_fullDict

    with open(outputFile, "w") as outFileP:
        for pos, meanDensity in sorted(finalDict.items()):
            outFileP.write(f"{pos}\t{meanDensity}\n")
