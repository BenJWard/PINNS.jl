module PINNS

using BioSequences
using RCall
using BioBridgeR.APE

export RUN_PINNS

function prepare_r(sequences::Vector{DNASequence}, seed::Int)
    R"""
    library(ape)
    library(phangorn)
    library(seqinr)

    set.seed($seed)

    SEQUENCES <- as.phyDat($sequences)

    PHYLOGENY <- pratchet(SEQUENCES, trace = 0)

    BASES <- c('a', 'c', 'g', 't')

    ancestryEstimation <- function(seqs, tree){
        anc <- ancestral.pars(tree, data, "ACCTRAN")
        m <- anc[[getRoot(tree)]]
        winners <- apply(m, 1, function(x) sample(BASES, size = 1, prob = x))
        return(as.DNAbin(winners[attr(anc, "index")]))
    }

    ANCESTOR <- ancestryEstimation(SEQUENCES, PHYLOGENY)
    """
end

function RUN_PINNS(sequences::Vector{DNASequence}, seed::Int = rand(0:10000000))
    println("The seed has been set to: ", seed, "...")
    println("Computing parsimony phylogeny, and computing ancestral sequence...")
    prepare_r(sequences, seed)
end

end # Module PINNS.
