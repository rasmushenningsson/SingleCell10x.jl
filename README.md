# SingleCell10x.jl

A Julia package for loading 10x Single Cell RNA-seq and antibody count data.

It supports loading of 10x ".h5" files as well as text-based ".mtx" and ".tsv/.csv" files.
The count matrix can be transposed at load time and the data types of the matrix and annotations can be chosen conveniently.


## Documentation
The package should hopefully be straightforward to use.
See below for example usage. More detailed descriptions are in the docstrings of the exported functions:

* `read10x`
* `read10x_matrix`
* `read10x_features`
* `read10x_barcodes`



## Example usage
Below, we will load files from the `test/` directory, but you can try with your own files instead.
```julia
basepath = joinpath(pkgdir(SingleCell10x),"test/data/500_PBMC_3p_LT_Chromium_X_50genes/")
```


### Basic loading
Loading a sparse count matrix, feature annotations and barcode annotations from a ".h5" file:
```julia
using SingleCell10x
filepath = joinpath(basepath, "filtered_feature_bc_matrix.h5")

matrix,features,barcodes = read10x(filepath)
```

You can also load the matrix/features/barcodes separately:
```julia
matrix2 = read10x_matrix(filepath)
features2 = read10x_features(filepath)
barcodes2 = read10x_barcodes(filepath)
```

### Sink types
Per default, the matrix will be a `SparseMatrixCSC{Int,Int}`, features will be a `NamedTuple` where keys are column names and elements are column vectors, and barcodes a `Vector{String}`.
But you can also specify sink types to load directly into the desired type:
```julia
using DataFrames
matrix,features,barcodes = read10x(filepath, Matrix{Float64}, DataFrame, DataFrame)
```
or similarly
```julia
matrix2 = read10x_matrix(filepath, Matrix{Float64})
features2 = read10x_features(filepath, DataFrame)
barcodes2 = read10x_barcodes(filepath, DataFrame)
```


### MatrixMarket .mtx files
Instead of loading a ".h5", `SingleCell10x.jl` also support the text-based "matrix.mtx(.gz)" with companion files. By default, feature/barcode filenames are searched for in the same folder.
```julia
filepath_mtx = joinpath(basepath,"filtered_feature_bc_matrix/matrix.mtx.gz")

matrix,features,barcodes = read10x(filepath_mtx)
```

If you load features/barcodes separately, you can either specify the path to the "matrix.mtx(.gz)" file or the feature/barcode file directly.

```julia
filepath_features = joinpath(basepath,"filtered_feature_bc_matrix/features.tsv.gz")
filepath_features = joinpath(basepath,"filtered_feature_bc_matrix/barcodes.tsv.gz")

features2 = read10x_features(filepath_mtx)
features3 = read10x_features(filepath_features)

barcodes2 = read10x_barcodes(filepath_mtx)
barcodes3 = read10x_barcodes(filepath_barcodes)
```
