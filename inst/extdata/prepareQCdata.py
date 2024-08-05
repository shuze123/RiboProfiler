import pysam
from collections import defaultdict

##############################################################################################################
# for mapping to genome qc analysis
##############################################################################################################

def prepareQCdata(longestTransInfo, samFile, outFile, seqType):
    geneinfoDict = {}

    with open(longestTransInfo, "r") as geneinfo:
        for line in geneinfo:
            # split tags
            _, gene_name, _, _, chr, strand, cdsRg, exon, utr5, cds, utr3 = line.split()
            utr5, cds, utr3 = int(utr5), int(cds), int(utr3)
            exonLength = utr5 + cds + utr3
            cdsStart, cdsEnd = (utr5 + 1), (utr5 + cds - 2)
            
            # get exon positions
            exonPositions = []
            for rg in exon.split(","):
                posSt, posEnd = map(int, rg.split(":"))
                exonPositions.extend(range(posSt, posEnd + 1))
            
            # calculate geneinfoDict
            counts = list(range(1, exonLength + 1))
            if strand == "+":
                pass
            else:
                counts.reverse()
            
            posKeys = [f"{chr}|{pos}" for pos in exonPositions]

            for i in range(len(posKeys)):
                geneinfoDict[posKeys[i]] = (cdsStart, cdsEnd, counts[i])

    print("Transforming genomic positions into transcriptome positions has been done successfully.")

    def RiboQcAnalysis(inputFile, outputFile, seq_type):
        frame_dict = defaultdict(int)

        reader = pysam.AlignmentFile(inputFile, "r")

        for record in reader:
            if not record.is_unmapped:
                refname = record.reference_name
                align_pos = record.reference_start
                read_length = record.query_length
                flag = record.flag

                if seq_type == "singleEnd":
                    if flag == 16:
                        end3Pos = align_pos + read_length - 1
                        readKey = f"{refname}|{end3Pos}"
                    elif flag == 0:
                        readKey = f"{refname}|{align_pos}"
                    else:
                        print("There are other flags!")
                        continue
                elif seq_type == "pairedEnd":
                    if flag == 16:
                        readKey = f"{refname}|{align_pos}"
                    elif flag == 0:
                        end3Pos = align_pos + read_length - 1
                        readKey = f"{refname}|{end3Pos}"
                    else:
                        print("There are other flags!")
                        continue

                if readKey in geneinfoDict:
                    start_codon_pos, stop_codon_pos, transPos = geneinfoDict[readKey]
                    rel2st = transPos - start_codon_pos
                    rel2sp = transPos - stop_codon_pos

                    frame_st = abs(rel2st) % 3
                    frame_sp = abs(rel2sp) % 3

                    if flag == 16:
                        align_pos_center = transPos + (read_length // 2)
                    elif flag == 0:
                        align_pos_center = transPos - (read_length // 2)

                    if align_pos_center <= start_codon_pos:
                        ftype = 2 # 5'UTR
                    elif start_codon_pos < align_pos_center <= stop_codon_pos + 2:
                        ftype = 3 # CDS
                    else:
                        ftype = 1 # 3'UTR

                    key = f"{read_length}\t{frame_st}\t{rel2st}\t{frame_sp}\t{rel2sp}\t{ftype}"
                    frame_dict[key] += 1

        reader.close()

        with open(outputFile, "w") as outfile:
            for key, val in frame_dict.items():
                outfile.write(f"{key}\t{val}\n")

    print("Processing sam files...")
    samFiles = samFile.split(",")
    outFiles = outFile.split(",")
    for i in range(len(samFiles)):
        RiboQcAnalysis(samFiles[i], outFiles[i], seq_type=seqType)
        print(f"{samFiles[i]} has been processed.")

##############################################################################################################
# for mapping to transcriptome qc analysis
##############################################################################################################

def prepareQCdata_ontrans(samFile, outFile, seqType):
    frame_dict = defaultdict(int)

    reader = pysam.AlignmentFile(samFile, "r")

    for record in reader:
        if not record.is_unmapped:
            refname = record.reference_name
            align_pos = record.reference_start
            read_length = record.query_length
            flag = record.flag

            if seqType == "singleEnd":
                if flag == 0:
                    exact_pos = align_pos
                elif flag == 16:
                    exact_pos = align_pos + read_length - 1
                else:
                    print("There are other flags!")
                    continue
            elif seqType == "pairedEnd":
                if flag == 0:
                    exact_pos = align_pos + read_length - 1
                elif flag == 16:
                    exact_pos = align_pos
                else:
                    print("There are other flags!")
                    continue

            start_codon_pos = int(refname.split("|")[4])
            stop_codon_pos = int(refname.split("|")[5])

            rel2st = exact_pos - start_codon_pos
            rel2sp = exact_pos - (stop_codon_pos + 1)

            frame_st = abs(rel2st) % 3
            frame_sp = abs(rel2sp) % 3

            if flag == 16:
                align_pos_center = exact_pos + (read_length // 2)
            elif flag == 0:
                align_pos_center = exact_pos - (read_length // 2)

            if align_pos_center < start_codon_pos:
                ftype = 2 # 5'UTR
            elif start_codon_pos <= align_pos_center <= stop_codon_pos:
                ftype = 3 # CDS
            else:
                ftype = 1 # 3'UTR

            key = f"{read_length}\t{frame_st}\t{rel2st}\t{frame_sp}\t{rel2sp}\t{ftype}"
            frame_dict[key] += 1

    reader.close()

    with open(outFile, "w") as outfile:
        for key, val in frame_dict.items():
            outfile.write(f"{key}\t{val}\n")
