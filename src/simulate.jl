

function nullsim(file::String, seed::Int, reps::Int, model::String, treemethod::String, ancmethod::String)
    sim_file = "nullsim_$file"
    anc_file = "ancestor_$file"
    prepare_r(seed)
    R"""
    sequences <- as.phyDat(read.dna($file, format = "fasta"))
    tree <- initialTree(sequences, $model, $treemethod)
    mlFit <- refineML(tree, sequences, $model)
    anc <- reconstructAncestor(mlFit, tree, $ancmethod)
    sims <- generateSimulations(mlFit, anc, $reps)
    rownames(sims) <- paste(names(sequences), rep(1:$reps, each = length(sequences)))
    sendToFile(sims, $sname)
    sendToFile(as.list(as.DNAbin(ancestor)), $aname)
    """
end

function prepare_r(seed::Int)
    R"""
    library(ape)
    library(phangorn)
    set.seed($seed)

    BASES <- c('a', 'c', 'g', 't')

    initialTree <- function(data, model, method){
        distance <- dist.ml(data, model)
        tree <- switch(tolower(method),
                       nj = NJ(distance),
                       upgma = upgma(distance),
                       stop("Invalid tree building method")
                )
        return(tree)
    }

    refineML <- function(tree, data, model){
        mlFit <- pml(tree, data, model)
        mlFit <- optim.pml(mlFit, model = model, control = pml.control(trace = 0))
        return(mlFit)
    }

    reconstructAncestor <- function(fitmodel, tree, method){
        anc <- ancestral.pml(fitmodel, tolower(method))
        ancMat <- anc[[getRoot(tree)]]
        winners <- apply(ancMat, 1, function(x) sample(BASES, size = 1, prob = x))
        ancestor <- winners[attr(anc, "index")]
        return(ancestor)
    }

    generateAlignment <- function(mlfit, root){
        return(as.DNAbin(simSeq(mlfit, l = length(root), rootseq = root)))
    }

    generateAlignments <- function(mlfit, root, reps){
        return(do.call(rbind, lapply(1:reps, function(x) generateSimulation(mlfit, root))))
    }

    sendToFile <- function(seqs, filename){
        write.dna(seqs, file = filename, format = "fasta")
    }
    """
end
