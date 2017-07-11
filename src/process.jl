
function process(realname::String, simname::String)
    output = open("$(basename(realname)).csv", "w")
    simreader = FASTA.Reader(open(simname, "r"))
    try
        println(output, "Seq1, Seq2, dN, dS, Real")
        sequence_records = open(realname, "r") do file
            collect(FASTA.Reader(file))
        end


        println(sequence_records)
        exit()
        write_rep_to_file(output, names, pairwise_dNdS(NG86, seqs), true)
        while !eof(simreader)
            for i in 1:endof(seqs)
        end
    finally
        close(output)
        close(simreader)
    end
end

@inline function write_rep_to_file(io, names, results, real)
    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        dN, dS = results[i, j]
        println(io, names[i], ", ", names[j], ", ", dN, ", ", dS, ", ", real)
    end
end

function readreal(filename)
    open(filename, "r") do file
        names = Vector{String}()
        sequences = Vector{DNASequence}()
        record = FASTA.Record()
        reader = FASTA.Reader(file)
        while !eof(reader)
            read!(reader, record)
            push!(names, FASTA.identifier(record))
            push!(sequences, FASTA.sequence(record))
        end
        return names, sequences
    end
end

function readsim(fh::FASTA.Reader)

end
