module Gating

using SingleCellProjections
using DataFrames
using Observables
using GLMakie
using UMAP

export gate

abstract type AbstractGater end

include("singlecell.jl")

end
