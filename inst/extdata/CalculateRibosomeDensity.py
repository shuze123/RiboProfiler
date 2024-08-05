import pysam

##############################################################################################################
# for mapping to genome density analysis
##############################################################################################################

def CalculateRibosomeDensity(inputFile, outputFile, min_length=23, max_length=35):
    # save in dictionary
    density_dict = {}
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
            read_length = record.query_length

            # filter read length
            if min_length <= read_length <= max_length:
                end5 = align_pos
                end3 = end5 + read_length - 1
                # shift +- 11nt
                centerEnd5 = end5 + 11
                centerEnd3 = end3 - 11
                centerLength = centerEnd3 - centerEnd5 + 1

                # ribo density
                for elem in range(centerEnd5, centerEnd3 + 1):
                    key = f"{refname}:{elem}"
                    density_dict[key] = density_dict.get(key, 0.0) + (1.0 / centerLength)

    reader.close()

    # output file
    with open(outputFile, "w") as outfile:
        # sort keys
        for key in sorted(density_dict.keys()):
            id, align_pos = key.split(":")
            raw_density = density_dict[key]

            # RPM normalization
            rpm = (raw_density / total_reads) * 1000000
            outfile.write(f"{id}\t{align_pos}\t{raw_density}\t{rpm}\n")

##############################################################################################################
# for mapping to transcriptome density analysis
##############################################################################################################

def CalculateRibosomeDensity_ontrans(inputFile, outputFile, min_length=23, max_length=35):
    # save in dictionary
    density_dict = {}
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
            read_length = record.query_length

            # filter read length
            if min_length < read_length < max_length:
                end5 = align_pos
                end3 = end5 + read_length - 1
                # shift +- 11nt
                centerEnd5 = end5 + 11
                centerEnd3 = end3 - 11
                centerLength = centerEnd3 - centerEnd5 + 1

                # ribo density
                for elem in range(centerEnd5, centerEnd3 + 1):
                    key = f"{refname}:{elem}"
                    density_dict[key] = density_dict.get(key, 0.0) + (1.0 / centerLength)

    reader.close()

    # output file
    with open(outputFile, "w") as outfile:
        # sort keys
        for key in sorted(density_dict.keys()):
            id, align_pos = key.split(":")
            raw_density = density_dict[key]

            # RPM normalization
            rpm = (raw_density / total_reads) * 1000000
            outfile.write(f"{id}\t{align_pos}\t{raw_density}\t{rpm}\n")
