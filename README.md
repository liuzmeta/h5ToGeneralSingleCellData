# h5ToGeneralSingleCellData
a script to convert h5 file to general single cell data
## Usage

e.g. `python convert_h5_to_general.py GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5` well create the following content in current directory:

- GSM3489182_Donor_01_filtered_gene_bc_matrices_h5
  - barcodes.tsv.gz
  - features.tsv.gz
  - matrix.mtx.gz