import pysam

def CalculateCountTPM(inputFile, outputFile, inputType):
    #####################################################
    # 1. count gene numbers
    #####################################################

    # save in dictionary
    gene_dict = {}

    # open sam file
    reader = pysam.AlignmentFile(inputFile, "r")

    # loop through each record
    for record in reader:
        # do something
        if not record.is_unmapped:
            # tags
            refname = record.reference_name

            # tags
            gene_name, _, _, cdsST, cdsSP, gene_length = refname.split("|")
            
            if inputType == "ribo":
                # tags
                align_pos = record.reference_start
                read_length = record.query_length

                # read center position
                align_pos_center = align_pos + (read_length // 2)
                
                # cds length
                cdsST, cdsSP = int(cdsST), int(cdsSP)
                cdsLength = cdsSP - cdsST + 1

                key = f"{gene_name}|{cdsLength}"
                # count reads ribo on CDS region
                if cdsST <= align_pos_center <= cdsSP:
                    if key not in gene_dict:
                        gene_dict[key] = 1
                    else:
                        gene_dict[key] += 1
            elif inputType == "rna":
                # count reads rna on transcript region
                key = f"{gene_name}|{gene_length}"
                if key not in gene_dict:
                    gene_dict[key] = 1
                else:
                    gene_dict[key] += 1
            else:
                print("error!")
                break

    reader.close()
    
    #####################################################
    # 2. calculate TPM values
    #####################################################

    # total reads one sample

    tpm_dict = {}

    # get TPM values
    for key, val in gene_dict.items():
        gene_name, gene_length = key.split("|")
        gene_length = int(gene_length)

        # reads / geneLength(kb)
        rpk = val / (gene_length / 1000)
        tpm_dict[gene_name] = [val, rpk]

    # output
    with open(outputFile, "w") as tpm_output:
        total_rpk = sum(rpk for count, rpk in tpm_dict.values())

        for key, val in tpm_dict.items():
            tpm = (val[1] / total_rpk) * 1000000
            # save
            tpm_output.write(f"{key}\t{val[0]}\t{tpm}\n")
