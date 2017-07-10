
function process(realname::String, simname::String)
    names, seqs = readreal(realname)
    println(names)
    println(seqs)
end

function readreal(fname::String)
    names = Vector{String}()
    sequences = Vector{DNASequence}()

    reader = FASTA.Reader(open(fname, "r"))
    record = FASTA.Record()

    while !eof(reader)
        read!(reader, record)
        push!(names, FASTA.identifier(record))
        push!(sequences, FASTA.sequence(record))
    end

    close(reader)

    return names, sequences
end
