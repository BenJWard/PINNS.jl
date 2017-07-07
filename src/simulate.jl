


function nullsim(file::String, seed::Int, reps::Int)
    prepare_r(seed)
    make_ancestry_r(file)
    generate_simulations(file, reps)
end

function prepare_r(seed::Int)
    R"""
    library(ape)
    library(phangorn)

    set.seed($seed)

    BASES <- c('a', 'c', 'g', 't')

    generateSimulation <- function(model, root){
        return(as.DNAbin(simSeq(model, l = length(root), rootseq = root)))
    }

    generateSimulations <- function(model, root, reps){
        return(do.call(rbind, lapply(1:reps, function(x) generateSimulation(model, root))))
    }

    putToFile <- function(seqs, filename){
        write.dna(seqs, file = filename, format = "fasta")
    }
    """
end

function load_seqs_r(file::String)
    R"""
    sequences <- as.phyDat(read.dna($file, format = "fasta"))
    """
end

function load_seqs_r(seqs)
    R"""
    sequences <- as.phyDat($seqs)
    """
end

function make_ancestry_r(file::String)
    load_seqs_r(file)
    make_ancestry_r()
end

function make_ancestry_r(seqs)
    load_seqs_r(seqs)
    make_ancestry_r()
end

function make_ancestry_r()
    R"""
    phylogeny <- acctran(pratchet(sequences, trace = 0), sequences)
    mlFit <- pml(phylogeny, sequences)
    mlFit <- optim.pml(mlFit, model="F81", control = pml.control(trace=0))
    anc <- ancestral.pml(mlFit, "ml")
    ancMat <- anc[[getRoot(phylogeny)]]
    winners <- apply(ancMat, 1, function(x) sample(BASES, size = 1, prob = x))
    ancestor <- winners[attr(anc, "index")]
    """
end

function generate_simulations(file::String, reps::Int)
    sname = "nullsim_$file"
    aname = "ancestor_$file"
    R"""
    sims <- generateSimulations(mlFit, ancestor, $reps)
    rownames(sims) <- paste(names(sequences), rep(1:$reps, each = length(sequences)))
    putToFile(sims, $sname)
    putToFile(as.list(as.DNAbin(ancestor)), $aname)
    """
end
