def GetGeneSinglePosDensity(geneInfo, inputFile, outputFile):
    ################################################################
    # 1. load density file
    ################################################################
    # save in dictionary
    density_dict = {}
    chr_dict = {}

    # open density file
    with open(inputFile, "r") as density_file:
        for line in density_file:
            # split tags
            chr, align_pos, _, rpm = line.split()
            density_dict[f"{chr}|{align_pos}"] = float(rpm)
            if chr not in chr_dict:
                chr_dict[chr] = 0

    ################################################################
    # 2. transform into transcriptome coordinate
    ################################################################
    with open(outputFile, "w") as out_gene_density:
        with open(geneInfo, "r") as gene_info:
            for line in gene_info:
                # split tags
                _, gene_name, gene_id, trans_id, chr, strand, _, exon, utr5, cds, utr3 = line.split()
                
                # exon length
                utr5, cds, utr3 = int(utr5), int(cds), int(utr3)

                # exclude no 5'UTR genes and can't be divided exactly by 3 genes
                exon_length = utr5 + cds + utr3

                # cds start and end position
                cds_start, cds_end = (utr5 + 1), (utr5 + cds - 2)

                # initialize counters
                if strand == "+":
                    count = 0
                    val = 1
                else:
                    count = exon_length + 1
                    val = -1

                # fill density
                for rg in exon.split(","):
                    pos_st, pos_end = map(int, rg.split(":")[1].split("-"))
                    for chr_pos in range(pos_st, pos_end + 1):
                        count += val
                        if chr in chr_dict:
                            density = density_dict.get(f"{chr}|{chr_pos}", 0)
                            out_gene_density.write(f"{gene_name}\t{trans_id}\t{count}\t{density}\n")

    # close file (handled by 'with' context manager)
