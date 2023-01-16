# Data source
10x example file "neurons_900_filtered_gene_bc_matrices.tar.gz" downloaded from here:
https://www.10xgenomics.com/resources/datasets/1-k-brain-cells-from-an-e-18-mouse-2-standard-2-0-1

# Code for subsetting the data matrix
```julia
using MatrixMarket
using DataFrames
using DelimitedFiles

ind = [22350, 12261, 15088, 7790, 27805, 8108, 15473, 25996, 8734, 21850, 25253, 11360, 15182, 27562, 26674, 1297, 25226, 21177, 4291, 18534, 17735, 14615, 20831, 21580, 27933, 22469, 9895, 17666, 14584, 116]

fn = "path/to/neurons_900_filtered_gene_bc_matrices/filtered_gene_bc_matrices/mm10/matrix.mtx"
X = MatrixMarket.mmread(fn)

fn_out = "matrix.mtx"
MatrixMarket.mmwrite(fn_out, X[ind,:])


fn_genes = "path/to/neurons_900_filtered_gene_bc_matrices/filtered_gene_bc_matrices/mm10/genes.tsv"
genes = readdlm(fn_genes, String)

fn_genes_out = "genes.tsv"
writedlm(fn_genes_out, genes[ind,:])

# barcode file used without changes


fn_dense_out = "dense.tsv"
writedlm(fn_dense_out,convert(Matrix,X[ind,:]))

```

