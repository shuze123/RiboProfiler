# load package
using XAM

# define function
function CalculateRNACoverage(;inputFile,outputFile,type = "coverage")
    # save in dict
    coverage_dict = Dict{String,Float64}()

    # open sam file
    reader = open(SAM.Reader,inputFile)
    record = SAM.Record()

    # loop
    while !eof(reader)
        empty!(record)
        read!(reader, record)
        # do something
        if SAM.ismapped(record)
            # tags
            refname = SAM.refname(record)
            align_pos = SAM.position(record)
            
            if type == "coverage" # get read coverage
                read_length = SAM.seqlength(record)
                # read right position
                End5 = align_pos + read_length - 1
                # ribo density
                for elem in range(align_pos,End5)
                    key = join([refname,elem],"|")
                    if !haskey(coverage_dict,key)
                        coverage_dict[key] = 1
                    else
                        coverage_dict[key] += 1
                    end
                end
            elseif type == "counts" # get read counts
                key = join([refname,align_pos],"|")
                if !haskey(coverage_dict,key)
                    coverage_dict[key] = 1
                else
                    coverage_dict[key] += 1
                end
            else
                println("error!")
                break
            end
        end
    end
    
    # output file
    outfile = open(outputFile,"w")

    # total densitys
    total_density = sum(values(coverage_dict))

    # sort keys
    for key in sort(collect(keys(coverage_dict)))
        id,align_pos = split(key,"|")
        raw_denisty = coverage_dict[key]
        
        # RPM normalization
        rpm = (raw_denisty/total_density)*1000000
        write(outfile,"$id\t$align_pos\t$raw_denisty\t$rpm\n")
    end
    close(outfile)
end