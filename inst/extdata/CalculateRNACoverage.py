import pysam

def CalculateRNACoverage(inputFile, outputFile, type="coverage"):
    # save in dictionary
    coverage_dict = {}
    total_reads = 0

    # open sam file
    reader = pysam.AlignmentFile(inputFile, "r")

    # loop through each record
    for record in reader:
        # do something
        if not record.is_unmapped:
            total_reads += 1
            # tags
            refname = record.reference_name
            align_pos = record.reference_start
            
            if type == "coverage": # get read coverage
                read_length = record.query_length
                # read right position
                End5 = align_pos + read_length - 1
                # ribo density
                for elem in range(align_pos, End5 + 1):
                    key = f"{refname}|{elem}"
                    # if key not in coverage_dict:
                    #     coverage_dict[key] = 1
                    # else:
                    #     coverage_dict[key] += 1
                    coverage_dict[key] = coverage_dict.get(key, 0) + 1
            elif type == "counts": # get read counts
                key = f"{refname}|{align_pos}"
                # if key not in coverage_dict:
                #     coverage_dict[key] = 1
                # else:
                #     coverage_dict[key] += 1
                coverage_dict[key] = coverage_dict.get(key, 0) + 1
            else:
                print("error!")
                break

    reader.close()

    # output file
    with open(outputFile, "w") as outfile:
        # sort keys
        for key in sorted(coverage_dict.keys()):
            id, align_pos = key.split("|")
            raw_density = coverage_dict[key]
            
            # RPM normalization
            rpm = (raw_density / total_reads) * 1000000
            outfile.write(f"{id}\t{align_pos}\t{raw_density}\t{rpm}\n")
