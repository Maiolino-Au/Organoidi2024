import argparse
import anndata
import scib
import scvi
import scanpy as sc
import muon as mu

parser = argparse.ArgumentParser(description="Convert and process H5AD file")
parser.add_argument("input_file", help="Path to the input file")
arg = parser.parse_args()

adata = sc.read_h5ad(arg.input_file)

adata.write(arg.input_file.replace(".h5ad", "_processed.h5ad"))