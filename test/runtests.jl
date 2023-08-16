using SingleCell10x
using SingleCell10x: RawIJV
using Test

using HDF5
using CodecZlib
using DelimitedFiles
using SparseArrays
using DataFrames


isgz(fn) = lowercase(splitext(fn)[2])==".gz"
_open(f, fn) = open(fn) do io
    f(isgz(fn) ? GzipDecompressorStream(io) : io)
end

read_matrix(fn,delim=',') = _open(io->readdlm(io,delim,Int), fn)
read_strings(fn,delim=',') = _open(io->readdlm(io,delim,String), fn)

function gunzip(source,dest)
    data = open(source,"r") do src
        read(GzipDecompressorStream(src))
    end
    write(dest, data)
end


@testset "SingleCell10x.jl" begin
    @testset "500_PBMC_50genes" begin
        base = "data/500_PBMC_3p_LT_Chromium_X_50genes/"
        filenames = ("filtered_feature_bc_matrix.h5", "filtered_feature_bc_matrix/matrix.mtx.gz")

        expected_mat = read_matrix(joinpath(base,"expected_matrix.csv"))
        expected_nnz = count(!iszero, expected_mat)
        expected_feature_ids = vec(read_strings(joinpath(base,"expected_feature_ids.csv")))
        expected_barcodes = vec(read_strings(joinpath(base,"expected_barcodes.csv")))

        expected_feature_names = read_strings(joinpath(base,"filtered_feature_bc_matrix/features.tsv.gz"),'\t')[:,2]
        expected_feature_types = fill("Gene Expression", 50)
        expected_feature_genome = fill("GRCh38", 50)
        expected_features_mtx = (;id=expected_feature_ids, name=expected_feature_names, feature_type=expected_feature_types)
        expected_features_h5 = (;expected_features_mtx..., genome=expected_feature_genome)

        m_fn = joinpath(base,"filtered_feature_bc_matrix/matrix.mtx.gz")
        f_fn = joinpath(base,"filtered_feature_bc_matrix/features.tsv.gz")
        b_fn = joinpath(base,"filtered_feature_bc_matrix/barcodes.tsv.gz")

        mtx_raw = split(open(x->readchomp(GzipDecompressorStream(x)), m_fn), '\n')

        @testset "$case" for (case,filename) in zip(("h5","mtx"), filenames)
            fn = joinpath(base,filename)

            expected_features = case=="mtx" ? expected_features_mtx : expected_features_h5

            X,f,b = read10x(fn)
            @test X == expected_mat
            @test X isa SparseMatrixCSC{Int,Int}
            @test f == expected_features
            @test b == expected_barcodes

            X,f,b = read10x(fn; transpose=true)
            @test X == expected_mat'
            @test X isa SparseMatrixCSC{Int,Int}
            @test f == expected_features
            @test b == expected_barcodes

            @test read10x_matrix_metadata(fn) == (size(expected_mat)..., expected_nnz)
            @test read10x_matrix_metadata(fn; transpose=true) == (size(expected_mat')..., expected_nnz)

            @test read10x_matrix(fn) == expected_mat
            @test read10x_matrix(fn; transpose=true) == expected_mat'

            @test read10x_features(fn) == expected_features
            @test read10x_barcodes(fn) == expected_barcodes
        end

        @testset "h5_extra" begin
            h5open(joinpath(base,filenames[1])) do h5
                X,f,b = read10x(h5)
                @test X == expected_mat
                @test X isa SparseMatrixCSC{Int,Int}
                @test f == expected_features_h5
                @test b == expected_barcodes

                X,f,b = read10x(h5; transpose=true)
                @test X == expected_mat'
                @test X isa SparseMatrixCSC{Int,Int}
                @test f == expected_features_h5
                @test b == expected_barcodes

                @test read10x_matrix_metadata(h5) == (size(expected_mat)..., expected_nnz)
                @test read10x_matrix_metadata(h5; transpose=true) == (size(expected_mat')..., expected_nnz)

                @test read10x_matrix(h5) == expected_mat
                @test read10x_matrix(h5; transpose=true) == expected_mat'

                @test read10x_features(h5) == expected_features_h5
                @test read10x_barcodes(h5) == expected_barcodes
            end
        end

        @testset "mtx_extra" begin
            @test read10x_features(m_fn) == expected_features_mtx
            @test_throws AssertionError read10x_features(m_fn; guess=nothing)
            @test read10x_features(f_fn) == expected_features_mtx
            @test read10x_features(f_fn; guess=nothing).id == expected_feature_ids

            @test read10x_barcodes(m_fn) == expected_barcodes
            @test read10x_barcodes(m_fn; guess=nothing) == mtx_raw # garbage, but we want to test that guess=nothing works
            @test read10x_barcodes(b_fn) == expected_barcodes
            @test read10x_barcodes(b_fn; guess=nothing) == expected_barcodes
        end

        expected_rawijv = RawIJV(findnz(sparse(expected_mat))..., size(expected_mat)...)
        expected_rawijv_t = RawIJV(expected_rawijv.J, expected_rawijv.I, expected_rawijv.V, expected_rawijv.N, expected_rawijv.P)

        matrix_sinks = (RawIJV, RawIJV{Float64}, RawIJV{Float64,Int}, SparseMatrixCSC, SparseMatrixCSC{Int32}, SparseMatrixCSC{Float64}, SparseMatrixCSC{Float64,Int32}, Matrix, Matrix{Int32}, Matrix{Float64})
        @testset "matrix_$(case)_$sink" for (case,filename) in zip(("h5","mtx"), filenames), sink in matrix_sinks
            fn = joinpath(base,filename)

            X = read10x_matrix(fn, sink)
            @test X isa sink
            if X isa RawIJV
                @test X==expected_rawijv
            else
                @test X==expected_mat
            end

            X2,f,b = read10x(fn, sink)
            @test X==X2

            X = read10x_matrix(fn, sink; transpose=true)
            @test X isa sink
            if X isa RawIJV
                @test X==expected_rawijv_t
            else
                @test X==expected_mat'
            end

            X2,f,b = read10x(fn, sink; transpose=true)
            @test X==X2
        end

        square(x) = x.^2
        barcode_sinks = (Vector, NamedTuple, DataFrame, square)
        barcode_ans = (expected_barcodes, (;barcode=expected_barcodes), DataFrame(barcode=expected_barcodes), expected_barcodes.^2)
        @testset "barcodes_$(case)_$sink" for (case,filename) in zip(("h5","mtx"), filenames), (sink,ans) in zip(barcode_sinks,barcode_ans)
            fn = joinpath(base,filename)

            b = read10x_barcodes(fn, sink)
            @test b isa (sink isa Function ? Vector : sink)
            @test b==ans

            X,f,b = read10x(fn, RawIJV, NamedTuple, sink)
            @test b isa (sink isa Function ? Vector : sink)
            @test b==ans
        end


        feature_sinks = (NamedTuple, DataFrame, identity)
        expected_ans = (expected_features_h5, expected_features_mtx)
        @testset "features_$(case)_$sink" for (case,filename,ans) in zip(("h5","mtx"), filenames, expected_ans), sink in feature_sinks
            fn = joinpath(base,filename)

            ans2 = sink==DataFrame ? DataFrame(ans) : ans

            f = read10x_features(fn, sink)
            @test f isa (sink isa Function ? NamedTuple : sink)
            @test f==ans2

            X,f,b = read10x(fn, RawIJV, sink)
            @test f isa (sink isa Function ? NamedTuple : sink)
            @test f==ans2
        end

        @testset "$(prefix)gunzipped" for prefix in ("","prefix_")
            mktempdir() do dir
                m_unzipped_fn = joinpath(dir, string(prefix,basename(splitext(m_fn)[1])))
                f_unzipped_fn = joinpath(dir, string(prefix,basename(splitext(f_fn)[1])))
                b_unzipped_fn = joinpath(dir, string(prefix,basename(splitext(b_fn)[1])))
                gunzip(m_fn, m_unzipped_fn)

                @test read10x_matrix(m_unzipped_fn) == expected_mat
                @test read10x_matrix(m_unzipped_fn; transpose=true) == expected_mat'

                @test_throws ErrorException read10x(m_unzipped_fn)
                @test_throws ErrorException read10x_features(m_unzipped_fn)
                @test_throws ErrorException read10x_barcodes(m_unzipped_fn)

                @test_throws AssertionError read10x_features(m_unzipped_fn; guess=nothing)
                @test read10x_barcodes(m_unzipped_fn; guess=nothing) == mtx_raw # garbage, but we want to test that guess=nothing works

                gunzip(f_fn, f_unzipped_fn)
                gunzip(b_fn, b_unzipped_fn)

                X,f,b = read10x(m_unzipped_fn)
                @test X==expected_mat
                @test f==expected_features_mtx
                @test b==expected_barcodes

                X,f,b = read10x(m_unzipped_fn; transpose=true)
                @test X==expected_mat'
                @test f==expected_features_mtx
                @test b==expected_barcodes

                @test read10x_features(m_unzipped_fn) == expected_features_mtx
                @test read10x_barcodes(m_unzipped_fn) == expected_barcodes

                @test read10x_features(f_unzipped_fn; guess=nothing) == expected_features_mtx
                @test read10x_barcodes(b_unzipped_fn; guess=nothing) == expected_barcodes
            end
        end

    end

    @testset "neurons_900" begin # Data from Cellranger 2.0.1
        base = "data/neurons_900_30genes/"
        m_fn = joinpath(base,"matrix.mtx")
        g_fn = joinpath(base,"genes.tsv")
        b_fn = joinpath(base,"barcodes.tsv")

        expected_mat = readdlm(joinpath(base,"dense.tsv"), '\t', Int)
        expected_features_raw = readdlm(joinpath(base,"expected_genes.tsv"), '\t', String)
        expected_features = (;id=expected_features_raw[:,1], name=expected_features_raw[:,2], feature_type=expected_features_raw[:,3])
        expected_barcodes = vec(readdlm(b_fn, '\t', String))

        mtx_raw = split(readchomp(m_fn), '\n')

        X,f,b = read10x(m_fn)
        @test X == expected_mat
        @test f == expected_features
        @test b == expected_barcodes

        X,f,b = read10x(m_fn; transpose=true)
        @test X == expected_mat'
        @test f == expected_features
        @test b == expected_barcodes

        @test read10x_matrix(m_fn) == expected_mat
        @test read10x_matrix(m_fn; transpose=true) == expected_mat'

        @test read10x_features(m_fn) == expected_features
        @test_throws AssertionError read10x_features(m_fn; guess=nothing)
        @test read10x_features(g_fn) == expected_features
        @test read10x_features(g_fn; guess=nothing) == expected_features

        @test read10x_barcodes(m_fn) == expected_barcodes
        @test read10x_barcodes(m_fn; guess=nothing) == mtx_raw # garbage, but we want to test that guess=nothing works
        @test read10x_barcodes(b_fn) == expected_barcodes
        @test read10x_barcodes(b_fn; guess=nothing) == expected_barcodes
    end
end
