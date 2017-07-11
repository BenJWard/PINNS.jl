
function process(realname::String, simname::String)
    names, seqs = readreal(realname)
    nseqs = length(seqs)
    println(names)
    println(seqs)
    output = open("$(basename(realname)).csv", "w")
    println(output, "Seq1, Seq2, dN, dS, Real")
    #simreader = FASTA.Reader(open(simname, "r"))
    write_rep_to_file(output, names, map(x -> x[1] / x[2], pairwise_dNdS(NG86, seqs)), true)


end

@inline function write_rep_to_file(io, names, results, real)
    @inbounds for i ∈ 1:endof(names), j ∈ (i + 1):endof(names)
        (dN, dS) = results[i, j]
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
