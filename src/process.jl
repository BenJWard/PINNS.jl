
function process(realname::String, simname::String)
    names, seqs = readreal(realname)
    nseqs = length(seqs)
    println(names)
    println(seqs)
    compute_dnds(seqs)
    simreader = FASTA.Reader(open(simname, "r"))
    output = open("$(basename(realname)).csv", "w")






end

function compute_dnds(seqs::Vector{DNASequence})
    for i ∈ 1:endof(seqs), j ∈ i:endof(seqs)
        icdns, jcdns = aligned_codons(seqs[i], seqs[j])
        println(icdns)
        println(jcdns)
    end
end

function aligned_codons{T<:NucleicAcidAlphabets}(x::BioSequence{T}, y::BioSequence{T}, start::Int = 1)
    xcdns = Vector{Kmer{eltype(T), 3}}()
    ycdns = Vector{Kmer{eltype(T), 3}}()
    pos = start
    while pos + 2 ≤ min(endof(x), endof(y))
        cdnx, okx = BioSequences.extract_kmer_impl(x, pos, 3)
        cdny, oky = BioSequences.extract_kmer_impl(y, pos, 3)
        if okx && oky
            push!(xcdns, convert(Kmer{eltype(T), 3}, cdnx))
            push!(ycdns, convert(Kmer{eltype(T), 3}, cdny))
        end
        pos += 3
    end
    return xcdns, ycdns
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
