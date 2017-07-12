
function process(realname::String, simname::String)
    output = open("$(basename(realname)).csv", "w")
    simreader = FASTA.Reader(open(simname, "r"))
    try
        println(output, "seq1, seq2, dN, dS, nmutations, real")
        sequence_records = open(realname, "r") do file
            collect(FASTA.Reader(file))
        end
        println("Processing real dataset...")
        process_seq_records(output, sequence_records, true)
        println("Processing simulated dataset...")
        while !eof(simreader)
            @inbounds for i in 1:endof(sequence_records)
                read!(simreader, sequence_records[i])
            end
            process_seq_records(output, sequence_records, false)
        end
        println("Finished.")
    finally
        close(output)
        close(simreader)
    end
end

function extract_first(d::Matrix{Tuple{Int,Int}})
    o = similar(d, Int)
    @inbounds for i in eachindex(d)
        o[i] = first(d[i])
    end
    return o
end

@inline function process_seq_records(output, records, isreal)
    sequences = FASTA.sequence.(records)
    mcounts = count_pairwise(Mutated, sequences...)
    write_rep_to_file(output,
                      FASTA.identifier.(records),
                      pairwise_dNdS(NG86, sequences),
                      mcounts,
                      isreal)
end

@inline function write_rep_to_file(io, names, results, counts, real)
    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        dN, dS = results[i, j]
        println(io, names[i], ", ", names[j], ", ", dN, ", ", dS, ", ", first(counts[i, j]), ", ", real)
    end
end
