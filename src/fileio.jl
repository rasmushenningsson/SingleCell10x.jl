struct AutoType end # Placeholder for using whatever eltype is stored in a file

struct RawIJV{Tv,Ti}
	I::Vector{Ti}
	J::Vector{Ti}
	V::Vector{Tv}
	P::Int
	N::Int
end

Base.:(==)(A::RawIJV, B::RawIJV) = A.P == B.P && A.N == B.N && A.I == B.I && A.J==B.J && A.V==B.V


LinearAlgebra.adjoint(X::RawIJV) = RawIJV(X.J, X.I, X.V, X.N, X.P)


Base.convert(::Type{SparseMatrixCSC{Tv,Ti}}, X::RawIJV{Tv,Ti}) where {Tv,Ti} = sparse(X.I, X.J, X.V, X.P, X.N)
Base.convert(::Type{SparseMatrixCSC{Tv}}, X::RawIJV{Tv,Ti}) where {Tv,Ti} = sparse(X.I, X.J, X.V, X.P, X.N)
Base.convert(::Type{SparseMatrixCSC}, X::RawIJV{Tv,Ti}) where {Tv,Ti} = sparse(X.I, X.J, X.V, X.P, X.N)

Base.convert(::Type{RawIJV{Tv,Ti}}, X::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = RawIJV(findnz(X)...,size(X)...)
Base.convert(::Type{RawIJV{Tv}}, X::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = RawIJV(findnz(X)...,size(X)...)
Base.convert(::Type{RawIJV}, X::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} = RawIJV(findnz(X)...,size(X)...)


Base.convert(::Type{Td}, X::Adjoint) where Td<:RawIJV = convert(Td, X.parent)'

function Base.convert(::Type{Matrix{T}}, X::RawIJV{T}) where T
	@assert length(X.I)==length(X.J)==length(X.V)
	Y = zeros(T, X.P, X.N)
	for (i,j,v) in zip(X.I, X.J, X.V)
		Y[i,j] += v
	end
	Y
end
Base.convert(::Type{Matrix}, X::RawIJV{T}) where T = convert(Matrix{T},X)



_ish5(fn) = lowercase(splitext(fn)[2]) == ".h5"



"""
	_open(f::Function, filename::AbstractString)

Opens file with .h5 and .gz autodetection.
"""
function _open(f::Function, filename::AbstractString)
	ext = lowercase(splitext(filename)[2])
	if ext==".h5"
		h5open(f, filename)
	else
		gz = lowercase(splitext(filename)[2])==".gz"
		open(filename) do io
			stream = gz ? GzipDecompressorStream(io) : io
			ret = f(stream)
			gz && close(stream)
			ret
		end
	end
end

function _delim(filename::AbstractString)
	fn,ext = splitext(lowercase(filename))
	ext == ".gz" && _delim(fn)
	ext == ".csv" ? ',' : '\t'
end



_matrix(::Type{T}, X, transpose::Bool) where T = convert(T, transpose ? X' : X)



function _read10x_matrix(io::HDF5.File, ::Type{Ti}, ::Type{Tv}) where {Ti,Tv}
	P,N,_ = _read10x_matrix_metadata(io)
	matrixGroup = HDF5.root(io)["matrix"]

	indptr = read(matrixGroup["indptr"])
	@assert length(indptr)==N+1
	rowval = read(matrixGroup["indices"])
	nzval = read(matrixGroup["data"])
	@assert length(rowval)==length(nzval)

	# convert (NB: this is a no-op if already of the correct type)
	if Ti !== AutoType
		indptr = convert(Vector{Ti}, indptr)
		rowval = convert(Vector{Ti}, rowval)
	end
	if Tv !== AutoType
		nzval = convert(Vector{Tv}, nzval)
	end

	rowval .+= 1 # 0-based to 1-based
	indptr .+= 1 # 0-based to 1-based

	SparseMatrixCSC(P, N, indptr, rowval, nzval)
end


function _read_mtx_lines!(I::Vector{Ti},J::Vector{Ti},V::Vector{Tv},nz,io) where {Ti,Tv}
	k = 0
	for line in eachline(io)
		k += 1
		k > nz && break

		i,j,v = eachsplit(line, ' ')
		I[k] = parse(Ti,i)
		J[k] = parse(Ti,j)
		V[k] = parse(Tv,v)
	end
	@assert k==nz "Mtx file doesn't have enough lines"
end

function _read_mtx_header(io)
	line = readline(io)
	m = match(r"^%%MatrixMarket matrix coordinate (real|integer) general$", line)
	@assert m !== nothing
	Tv = m.captures[1]=="real" ? Float64 : Int

	# skip comments
	while isempty(line) || line[1]=='%'
		line = readline(io)
	end
	P,N,nz = parse.(Int, split(line))
	Tv,P,N,nz
end


function _read10x_matrix(io, ::Type{Ti}, ::Type{Tv}) where {Ti,Tv}
	Tv_file,P,N,nz = _read_mtx_header(io)

	I = zeros(Ti !== AutoType ? Ti : Int, nz)
	J = zeros(Ti !== AutoType ? Ti : Int, nz)
	V = zeros(Tv !== AutoType ? Tv : Tv_file, nz)
	_read_mtx_lines!(I,J,V,nz,io)

	RawIJV(I,J,V,P,N)
end

_read10x_matrix(fn::AbstractString, ::Type{Ti}, ::Type{Tv}) where {Ti,Tv} = _open(x->_read10x_matrix(x,Ti,Tv), fn)



_read10x_matrix(::Type{SparseMatrixCSC{Tv,Ti}}, args...) where {Tv,Ti} = _read10x_matrix(args..., Ti, Tv)
_read10x_matrix(::Type{SparseMatrixCSC{Tv}}, args...) where Tv = _read10x_matrix(args..., AutoType, Tv)
_read10x_matrix(::Type{SparseMatrixCSC}, args...) = _read10x_matrix(args..., AutoType, AutoType)

_read10x_matrix(::Type{RawIJV{Tv,Ti}}, args...) where {Tv,Ti} = _read10x_matrix(args..., Ti, Tv)
_read10x_matrix(::Type{RawIJV{Tv}}, args...) where Tv = _read10x_matrix(args..., AutoType, Tv)
_read10x_matrix(::Type{RawIJV}, args...) = _read10x_matrix(args..., AutoType, AutoType)

_read10x_matrix(::Type{Matrix{Tv}}, args...) where Tv = _read10x_matrix(args..., AutoType, Tv)
_read10x_matrix(::Type{Matrix}, args...) = _read10x_matrix(args..., AutoType, AutoType)



"""
	read10x_matrix(io, [matrixtype=SparseMatrixCSC{Int,Int}; transpose=false)

Read 10x (CellRanger) count matrix.
`io` can be filename (".h5" or ".mtx[.gz]") or an `HDF5.File` handle or an IO object.

The returned count matrix is of type `matrixtype` and can be one of:
* `SparseMatrixCSC` - The standard sparse matrix format in Julia, found in `SparseArrays`.
* `Matrix` - A standard dense matrix.

If `transpose` is `true`, the data matrix will be read transposed.

See also: [`read10x`](@ref), [`read10x_barcodes`](@ref), [`read10x_features`](@ref)
"""
function read10x_matrix(io, ::Type{T}=SparseMatrixCSC{Int,Int}; transpose=false) where T
	X = _read10x_matrix(T, io)
	_matrix(T,X,transpose)
end




function _read10x_matrix_metadata(io::HDF5.File)
	matrixGroup = HDF5.root(io)["matrix"]
	shape = read(matrixGroup["shape"])
	@assert length(shape)==2
	P = convert(Int,shape[1])
	N = convert(Int,shape[2])
	data = matrixGroup["data"]
	nz = length(data)
	P,N,nz
end

function _read10x_matrix_metadata(io)
	_,P,N,nz = _read_mtx_header(io)
	P,N,nz
end

_read10x_matrix_metadata(fn::AbstractString) = _open(_read10x_matrix_metadata, fn)

"""
	read10x_matrix_metadata(io; transpose=false)

Read 10x (CellRanger) count matrix metadata.
`io` can be filename (".h5" or ".mtx[.gz]") or an `HDF5.File` handle or an IO object.

Returns `(P,N,nnz)` the height, width and number of nonzero entries of the count matrix.

If `transpose` is `true`, then `P` and `N` are swapped.

See also: [`read10x`](@ref), [`read10x_matrix`](@ref)
"""
function read10x_matrix_metadata(io; transpose=false)::Tuple{Int,Int,Int}
	P,N,nz = _read10x_matrix_metadata(io)
	transpose ? (N,P,nz) : (P,N,nz)
end



function _guessfilename(filename::AbstractString, names::Tuple, description)
	dir,base = splitdir(filename)

	pattern = r"matrix\.mtx(\.gz)?$"
	occursin(pattern, base) || return filename

	for name in names, ext in (".tsv.gz", ".tsv", ".csv.gz", ".csv")
		fn = joinpath(dir, replace(base, pattern=>string(name,ext)))
		isfile(fn) && return fn
	end

	prefix = length(names)==1 ? names[1] : string('(',join(names,'|'),')')
	error("Could not find $description file matching \"$filename\". Tried $prefix.(tsv|csv)[.gz].")
end


guessfeaturefilename(fn) = _guessfilename(fn, ("features","genes"), "feature")
guessbarcodefilename(fn) = _guessfilename(fn, ("barcodes",), "barcode")



function _read10x_barcodes(io::HDF5.File)
	read(HDF5.root(io)["matrix"]["barcodes"])
end
function _read10x_barcodes(io; delim)
	barcodes = readdlm(io, delim, String)
	@assert size(barcodes,2)==1
	barcode=vec(barcodes)
end
_read10x_barcodes(fn::AbstractString; kwargs...) = _open(fn) do io
	_read10x_barcodes(io; kwargs...)
end


_read10x_barcodes_autodetect(io;kwargs...) = _read10x_barcodes(io;kwargs...)
function _read10x_barcodes_autodetect(filename::AbstractString; kwargs...)
	if _ish5(filename)
		_read10x_barcodes(filename; kwargs...)
	else
		_read10x_barcodes_triplet(filename; kwargs...)
	end
end
function _read10x_barcodes_triplet(filename; guess=guessbarcodefilename, kwargs...)
	fn = _filename(guess,filename)
	_read10x_barcodes(fn; delim=_delim(fn), kwargs...)
end



_barcodes(::Type{T}, b) where T = T((;barcode=b))
_barcodes(::Type{T}, b) where {T<:AbstractVector} = convert(T,b)
_barcodes(fun, b) = fun(b)



"""
	read10x_barcodes(io, [barcodetype=Vector]; [guess, delim])

Read 10x barcodes.
`io` can be a filename (".h5", ".tsv[.gz]", ".csv[.gz]" or "[prefix]matrix.mtx[.gz]") or an `HDF5.File` handle or an IO object.

The returned annotation object for barcodes is of the type `barcodetype` and can be one of:
* Vector - With barcode values for each cell.
* NamedTuple - With "barcode" as key, and barcode vector as value.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple as input.

If the input filename is "[prefix]matrix.mtx[.gz]", then the same folder will be searched for a matching barcode file.
Set `guess=nothing` to disable. `guess` can also be set to a function taking the matrix file path as input and returning the corresponding barcode file path.

The column delimiter `delim` will be detected from the file name if not specifed.


See also: [`read10x`](@ref), [`read10x_matrix`](@ref), [`read10x_features`](@ref)
"""
function read10x_barcodes(io, barcodetype=Vector; kwargs...)
	_barcodes(barcodetype, _read10x_barcodes_autodetect(io;kwargs...))
end



function _read10x_features(io::HDF5.File)
	featureGroup = HDF5.root(io)["matrix"]["features"]
	cols = vcat("id", "name", "feature_type", read(featureGroup["_all_tag_keys"]))
	(; map(x->(Symbol(x), read(featureGroup,x)), cols)...)
end
function _read10x_features(io; delim='\t')
	featuresRaw = readdlm(io, delim, String)
	@assert size(featuresRaw,2)>=2

	id = featuresRaw[:,1]
	name = featuresRaw[:,2]
	feature_type=size(featuresRaw,2)>=3 ? featuresRaw[:,3] : fill("Gene Expression",length(id))
	(;id, name, feature_type)
end
_read10x_features(fn::AbstractString; kwargs...) = _open(fn) do io
	_read10x_features(io; kwargs...)
end


_read10x_features_autodetect(io;kwargs...) = _read10x_features(io;kwargs...)
function _read10x_features_autodetect(filename::AbstractString; kwargs...)
	if _ish5(filename)
		_read10x_features(filename; kwargs...)
	else
		_read10x_features_triplet(filename; kwargs...)
	end
end
function _read10x_features_triplet(filename; guess=guessfeaturefilename, kwargs...)
	fn = _filename(guess,filename)
	_read10x_features(fn; delim=_delim(fn), kwargs...)
end



_features(fun, f) = fun(f)


"""
	read10x_features(io, [featuretype=NamedTuple]; [guess, delim])

Read 10x features.
`io` can be a filename (".h5", ".tsv[.gz]", ".csv[.gz]" or ".mtx[.gz]") or an `HDF5.File` handle or an IO object.

The returned annotation table for features is of the type `featuretype` and can be one of:
* NamedTuple - With column names as keys and column vectors as values.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple as input.

If the input filename is "[prefix]matrix.mtx[.gz]", then the same folder will be searched for a matching feature file.
Set `guess=nothing` to disable. `guess` can also be set to a function taking the matrix file path as input and returning the corresponding feature file path.

The column delimiter `delim` will be detected from the file name if not specifed.

See also: [`read10x`](@ref), [`read10x_matrix`](@ref), [`read10x_barcodes`](@ref)
"""
function read10x_features(io, featuretype=NamedTuple; kwargs...)
	_features(featuretype, _read10x_features_autodetect(io;kwargs...))
end




function _read10x(::Type{T}, io::HDF5.File) where T
	f = _read10x_features(io)
	b = _read10x_barcodes(io)
	X = _read10x_matrix(T, io) # read matrix last to return quicker if features/barcodes error
	X,f,b
end
_read10x(::Type{T}, filename::AbstractString) where T = h5open(x->_read10x(T,x), filename)

function _read10x(::Type{T}, matrix_filename::AbstractString, features_filename::AbstractString, barcodes_filename::AbstractString; kwargs...) where T
	f = _read10x_features(features_filename; delim=_delim(features_filename), kwargs...)
	b = _read10x_barcodes(barcodes_filename; delim=_delim(barcodes_filename), kwargs...)
	X = _read10x_matrix(T, matrix_filename) # read matrix last to return quicker if features/barcodes error
	X,f,b
end



_filename(f::AbstractString, ::Any) = f
_filename(f, filename::AbstractString) = f(filename)
_filename(::Nothing, filename::AbstractString) = filename


_read10x_autodetect(::Type{T}, io::HDF5.File; kwargs...) where T = _read10x(T, io; kwargs...)
function _read10x_autodetect(::Type{T}, filename::AbstractString; kwargs...) where T
	if _ish5(filename)
		_read10x(T, filename; kwargs...)
	else
		_read10x_triplet(T, filename; kwargs...)
	end
end
function _read10x_triplet(::Type{T}, filename; features=guessfeaturefilename, barcodes=guessbarcodefilename, kwargs...) where T
	features = _filename(features, filename)
	barcodes = _filename(barcodes, filename)
	features == filename && error("Failed to guess feature file name")
	barcodes == filename && error("Failed to guess barcode file name")
	_read10x(T, filename, features, barcodes; kwargs...)
end



"""
	read10x(io, [matrixtype=SparseMatrixCSC{Int,Int}, featuretype=NamedTuple, barcodetype=Vector];
	        transpose=false, [features, barcodes, delim])

Read 10x (CellRanger) count matrix and annotations.
`io` can be filename (".h5" or "[prefix]matrix.mtx[.gz]") or an `HDF5.File` handle.
NB: If the input filename is "[prefix]matrix.mtx[.gz]", the same folder will be searched for matching feature and barcode files.

The returned count matrix is of type `matrixtype` and can be one of:
* `SparseMatrixCSC` - The standard sparse matrix format in Julia, found in `SparseArrays`.
* `Matrix` - A standard dense matrix.

The returned annotations for barcodes/features are of types `barcodetype`/`featuretype` and can be one of:
* NamedTuple - With column names as keys and column vectors as values.
* Vector - For barcodes only.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple (features) or Vector (barcodes) as input.

If `transpose` is `true`, the data matrix will be read transposed.

`features` and `barcodes` are only allowed for `mtx[.gz]` files.
They can be used to specify the file paths for `features`/`barcodes`.
By default, the file paths will be autodetected.
Set to `nothing` to disable.
They can also be set to functions taking the matrix file path as input and returning the corresponding barcode/feature file path.

The column delimiter `delim` is only used for features/barcodes.(tsv|csv) files, and will be detected from the file name if not specifed.

See also: [`read10x_matrix`](@ref), [`read10x_barcodes`](@ref), [`read10x_features`](@ref)
"""
function read10x(io, ::Type{T}=SparseMatrixCSC{Int,Int}, featuretype=NamedTuple, barcodetype=Vector;
                 transpose=false, kwargs...) where T
	X, f, b = _read10x_autodetect(T, io; kwargs...)

	P,N = X isa RawIJV ? (X.P,X.N) : size(X)

	i = findfirst(col->length(col)!=P, values(f))
	@assert i===nothing "Inconsistent number of features, expected $P rows, got $(length(f[i]))."
	N2 = length(b)
	@assert N2 == N "Inconsistent number of barcodes, expected $N rows, got $N2."

	_matrix(T,X,transpose), _features(featuretype,f), _barcodes(barcodetype,b)
end
