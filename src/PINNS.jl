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
    end
    return parse_args(s)
end

function cmd()
    arguments = parse_command_line()
    if arguments["%COMMAND%"] == "nullsim"
        args = arguments["nullsim"]
        nullsim(args["inputfile"], args["seed"], args["reps"])
    end
end

include("simulate.jl")

end # Module PINNS.
