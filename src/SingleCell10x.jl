module SingleCell10x

export
	read10x_matrix,
	read10x_matrix_metadata,
	read10x_features,
	read10x_barcodes,
	read10x

using DelimitedFiles
using HDF5
using LinearAlgebra
using SparseArrays
using CodecZlib

include("fileio.jl")

end
