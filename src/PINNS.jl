module PINNS

using ArgParse
using BioSequences
using RCall
using BioBridgeR.APE

export cmd

function parse_command_line()
    s = ArgParseSettings()
    s.add_help = true

    @add_arg_table s begin
        "nullsim"
            help = "Simulate null distribution."
            action = :command
    end

    @add_arg_table s["nullsim"] begin
        "inputfile"
            help = "An input fasta file."
            arg_type = String
            required = true
        "reps"
            help = "Number of simulations to produce."
            arg_type = Int64
            required = true
        "seed"
            help = "A seed to set for RNG"
            arg_type = Int64
            default = rand(1:1000000)
        "--treemethod", "-t"
            help = "Method used to construct the starting tree - 'nj' or 'upgma'."
            arg_type = String
            default = "nj"
        "--submodel", "-s"
            help = "The substitution model used in model fitting and ancestral sequence estimation - 'JC' or 'F81'."
            arg_type = String
            default = "JC"
        "--ancmethod", "-a"
            help = "the criterion to assign the internal nodes during ancestral sequence reconstruction. - 'marginal', 'ml', or 'bayes'."
            arg_type = String
            default = "ml"
    end
    return parse_args(s)
end

function cmd()
    arguments = parse_command_line()
    if arguments["%COMMAND%"] == "nullsim"
        args = arguments["nullsim"]
        nullsim(args["inputfile"], args["seed"], args["reps"], args["submodel"], args["treemethod"], args["ancmethod"])
    end
end

include("simulate.jl")

end # Module PINNS.
