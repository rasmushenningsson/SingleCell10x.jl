const RawCSC{Tv,Ti} = Tuple{Vector{Ti},Vector{Ti},Vector{Tv},Int,Int}
_transpose((I,J,V,P,N)::RawCSC) = (J,I,V,N,P)

function _ismatrixmtx(fn)
	fn = lowercase(basename(fn))
	fn=="matrix.mtx" || fn=="matrix.mtx.gz"
end


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

function _delimkwarg(filename::AbstractString)
	fn,ext = splitext(lowercase(filename))
	ext == ".h5" && return (;)
	if ext==".gz"
		fn,ext=splitext(fn)
	end
	ext == ".csv" && return (;delim=',')
	(;delim='\t')
end




function _matrix(::Type{RawCSC{Tv,Ti}}, (I,J,V,P,N)) where {Tv,Ti}
	convert(Vector{Ti},I), convert(Vector{Ti},J), convert(Vector{Tv},V), convert(Int,P), convert(Int,N)
end

_matrix(::Type{RawCSC{Tv}}, X) where Tv = _matrix(RawCSC{Tv,Int}, X)
_matrix(::Type{RawCSC}, X) = _matrix(RawCSC{Int,Int}, X)


_matrix(::Type{SparseMatrixCSC{Tv,Ti}}, X) where {Tv,Ti} = sparse(_matrix(RawCSC{Tv,Ti}, X)...)
_matrix(::Type{SparseMatrixCSC{Tv}}, X) where Tv = sparse(_matrix(RawCSC{Tv}, X)...)
_matrix(::Type{SparseMatrixCSC}, X) = sparse(_matrix(RawCSC, X)...)


function _matrix(::Type{Matrix{T}}, (I,J,V,P,N)::RawCSC{Tv}) where {T,Tv}
	@assert length(I)==length(J)==length(V)
	Y = zeros(T, P, N)
	for (i,j,v) in zip(I,J,V)
		Y[i,j] += v
	end
	Y
end
_matrix(::Type{Matrix}, X::RawCSC{Tv}) where Tv = _matrix(Matrix{Tv}, X)



# Fallback for user-defined sinks - f should be a function or constructor taking I,J,V,P,N
_matrix(f, X) = f(X...)


function _matrix(matrixtype, X, transpose::Bool)
	_matrix(matrixtype, transpose ? _transpose(X) : X)
end





function _read10x_matrix(io::HDF5.File)::RawCSC
	P,N,_ = _read10x_matrix_metadata(io)
	matrixGroup = HDF5.root(io)["matrix"]

	I = read(matrixGroup["indices"]) .+ 1 # 0-based to 1-based
	indptr = read(matrixGroup["indptr"])
	@assert length(indptr)==N+1
	J = zeros(Int,length(I))
	for k=1:N
		J[ indptr[k]+1:indptr[k+1] ] .= k
	end
	V = read(matrixGroup["data"])

	I,J,V,P,N
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


function _read10x_matrix(io)::RawCSC
	Tv,P,N,nz = _read_mtx_header(io)

	I = zeros(Int,nz)
	J = zeros(Int,nz)
	V = zeros(Tv,nz)
	_read_mtx_lines!(I,J,V,nz,io)

	I,J,V,P,N
end

_read10x_matrix(fn::AbstractString) = _open(_read10x_matrix, fn)



"""
	read10x_matrix(io, [matrixtype=SparseMatrixCSC{Int,Int}; transpose=false)

Read 10x (CellRanger) count matrix.
`io` can be filename (".h5" or ".mtx(.gz)") or an `HDF5.File` handle or an IO object.

The returned count matrix is of type `matrixtype` and can be one of:
* `SparseMatrixCSC` - The standard sparse matrix format in Julia, found in `SparseArrays`.
* `RawCSC` - A tuple of `I,J,V,P,N`, suitable for passing to the `sparse` function in `SparseArrays`.
* `Matrix` - A standard dense matrix.
* A user defined function/type constructor taking `I,J,V,P,N` as input.

If `transpose` is `true`, the data matrix will be read transposed.

See also: [`read10x`](@ref), [`read10x_barcodes`](@ref), [`read10x_features`](@ref)
"""
function read10x_matrix(io, matrixtype=SparseMatrixCSC{Int,Int}; transpose=false)
	X = _read10x_matrix(io)
	_matrix(matrixtype,X,transpose)
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
`io` can be filename (".h5" or ".mtx(.gz)") or an `HDF5.File` handle or an IO object.

Returns `(P,N,nnz)` the height, width and number of nonzero entries of the count matrix.

If `transpose` is `true`, then `P` and `N` are swapped.

See also: [`read10x`](@ref), [`read10x_matrix`](@ref)
"""
function read10x_matrix_metadata(io; transpose=false)::Tuple{Int,Int,Int}
	P,N,nz = _read10x_matrix_metadata(io)
	transpose ? (N,P,nz) : (P,N,nz)
end



function _guessfilename(filename::AbstractString, names::Tuple, description)
	_ismatrixmtx(filename) || return filename
	dir,base = splitdir(filename)

	for name in names, ext in (".tsv.gz", ".tsv", ".csv.gz", ".csv")
		fn = joinpath(dir, string(name,ext))
		isfile(fn) && return fn
	end

	prefix = length(names)==1 ? names[1] : string('(',join(names,'|'),')')
	error("Could not find $description file matching \"$filename\". Tried $prefix.(tsv|csv)[.gz].")
end




function _read10x_barcodes(io::HDF5.File)
	read(HDF5.root(io)["matrix"]["barcodes"])
end
function _read10x_barcodes(io; delim)
	barcodes = readdlm(io, delim, String)
	@assert size(barcodes,2)==1
	barcode=vec(barcodes)
end

function _read10x_barcodes(fn::AbstractString; guessfilename=true, kwargs...)
	if guessfilename
		fn = _guessfilename(fn, ("barcodes",), "barcode")
	end

	_open(fn) do io
		_read10x_barcodes(io; _delimkwarg(fn)..., kwargs...)
	end
end


_barcodes(::Type{T}, b) where T = T((;barcode=b))
_barcodes(::Type{T}, b) where {T<:AbstractVector} = convert(T,b)
_barcodes(fun, b) = fun(b)



"""
	read10x_barcodes(io, [barcodetype=Vector]; delim, guessfilename=true)

Read 10x barcodes.
`io` can be a filename (".h5", ".tsv(.gz)", ".csv(.gz)" or "matrix.mtx(.gz)") or an `HDF5.File` handle or an IO object.

The returned annotation object for barcodes is of the type `barcodetype` and can be one of:
* Vector - With barcode values for each cell.
* NamedTuple - With "barcode" as key, and barcode vector as value.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple as input.

The column delimiter `delim` will be guessed from the file name if not specifed.

If `guessfilename` is true and the input filename is "matrix.mtx(.gz)", then the same folder will be searched for a matching barcode file.

See also: [`read10x`](@ref), [`read10x_matrix`](@ref), [`read10x_features`](@ref)
"""
function read10x_barcodes(io, barcodetype=Vector; kwargs...)
	_barcodes(barcodetype, _read10x_barcodes(io;kwargs...))
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


function _read10x_features(fn::AbstractString; guessfilename=true, kwargs...)
	if guessfilename
		fn = _guessfilename(fn, ("features","genes"), "feature")
	end

	_open(fn) do io
		_read10x_features(io; _delimkwarg(fn)..., kwargs...)
	end
end


_features(fun, f) = fun(f)



"""
	read10x_features(io, [featuretype=NamedTuple]; delim, guessfilename=true)

Read 10x features.
`io` can be a filename (".h5", ".tsv(.gz)", ".csv(.gz)" or ".mtx(.gz)") or an `HDF5.File` handle or an IO object.

The returned annotation table for features is of the type `featuretype` and can be one of:
* NamedTuple - With column names as keys and column vectors as values.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple as input.

The column delimiter `delim` will be guessed from the file name if not specifed.

If `guessfilename` is true and the input filename is "matrix.mtx(.gz)", then the same folder will be searched for a matching feature file.

See also: [`read10x`](@ref), [`read10x_matrix`](@ref), [`read10x_barcodes`](@ref)
"""
function read10x_features(io, featuretype=NamedTuple; kwargs...)
	_features(featuretype, _read10x_features(io;kwargs...))
end




function _read10x(io::HDF5.File)
	f = _read10x_features(io)
	b = _read10x_barcodes(io)
	X = _read10x_matrix(io) # read matrix last to return quicker if features/barcodes error
	X,f,b
end
function _read10x(filename::AbstractString; kwargs...)
	if lowercase(splitext(filename)[2]) == ".h5"
		h5open(filename) do io # only open the file once
			_read10x(io; kwargs...)
		end
	else
		# TODO: check that filename is matrix.mtx.(.gz) because otherwise we cannot find features/barcodes
		_ismatrixmtx(filename) || error("""Expected filename to end with ".h5" or be "matrix.mtx(.gz)".""")
		f = _read10x_features(filename; kwargs...)
		b = _read10x_barcodes(filename; kwargs...)
		X = _read10x_matrix(filename) # read matrix last to return quicker if features/barcodes error
		X,f,b
	end
end



"""
	read10x(io, [matrixtype=SparseMatrixCSC{Int,Int}, featuretype=NamedTuple, barcodetype=Vector]; transpose=false, delim)

Read 10x (CellRanger) count matrix and annotations.
`io` can be filename (".h5" or "matrix.mtx(.gz)") or an `HDF5.File` handle.
NB: If the input filename is "matrix.mtx(.gz)", the same folder will be searched for matching feature and barcode files.

The returned count matrix is of type `matrixtype` and can be one of:
* `SparseMatrixCSC` - The standard sparse matrix format in Julia, found in `SparseArrays`.
* `RawCSC` - A tuple of `I,J,V,P,N`, suitable for passing to the `sparse` function in `SparseArrays`.
* `Matrix` - A standard dense matrix.
* A user defined function/type constructor taking `I,J,V,P,N` as input.

The returned annotations for barcodes/features are of types `barcodetype`/`featuretype` and can be one of:
* NamedTuple - With column names as keys and column vectors as values.
* Vector - For barcodes only.
* DataFrame - See `DataFrames` package.
* Other table types that can be constructed from NamedTuples.
* A user defined function/type constructor taking a NamedTuple as input.

If `transpose` is `true`, the data matrix will be read transposed.

The column delimiter `delim` is only used for features/barcodes.(tsv|csv) files, and will be guessed from the file name if not specifed.

See also: [`read10x_matrix`](@ref), [`read10x_barcodes`](@ref), [`read10x_features`](@ref)
"""
function read10x(io, matrixtype=SparseMatrixCSC{Int,Int}, featuretype=NamedTuple, barcodetype=Vector; transpose=false, kwargs...)
	(I,J,V,P,N), f, b = _read10x(io; kwargs...)

	i = findfirst(col->length(col)!=P, values(f))
	@assert i===nothing "Inconsistent number of features, expected $P rows, got $(length(f[i]))."
	N2 = length(b)
	@assert N2 == N "Inconsistent number of barcodes, expected $N rows, got $N2."

	_matrix(matrixtype,(I,J,V,P,N),transpose), _features(featuretype,f), _barcodes(barcodetype,b)
end
