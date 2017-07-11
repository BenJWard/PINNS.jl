
function process(realname::String, simname::String)
    output = open("$(basename(realname)).csv", "w")
    simreader = FASTA.Reader(open(simname, "r"))
    try
        println(output, "Seq1, Seq2, dN, dS, Real")
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

@inline function process_seq_records(output, records, isreal)
    write_rep_to_file(output,
                      FASTA.identifier.(records),
                      pairwise_dNdS(NG86,
                      FASTA.sequence.(records)),
                      isreal)
end

@inline function write_rep_to_file(io, names, results, real)
    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        dN, dS = results[i, j]
        println(io, names[i], ", ", names[j], ", ", dN, ", ", dS, ", ", real)
    end
end
