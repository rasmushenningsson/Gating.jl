module Gating

using SingleCellProjections
using DataFrames
using Observables
using GLMakie
using Colors
using UMAP
using PointInPoly

export gate,
       annotate!

abstract type AbstractGater end

include("singlecell.jl")

end
