module PINNS

using BioSequences
using RCall
using BioBridgeR.APE

function prepare_r_session()
    R"""
    library(ape)
    library(phangorn)
    """
end

end # Module PINNS.
