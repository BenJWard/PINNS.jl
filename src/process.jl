
function process(realname::String, simname::String)
    names, seqs = readreal(realname)
    nseqs = length(seqs)
    println(names)
    println(seqs)




    simreader = FASTA.Reader(open(simname), "r")




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
