import gzip
import os
import scipy.sparse as sp_sparse
import scipy.io
import tables
import sys
import io


def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'GRCh38')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        genes = f.get_node(mat_group, 'genes').read()
        gene_names = f.get_node(mat_group, 'gene_names').read()

        return (genes, gene_names, barcodes, matrix)


def writefile(fname, data):
    genes, gene_names, barcodes, matrix = data
    sample = '_'.join(fname.split('.')[:-1])
    os.mkdir(sample)
    b_s = io.BytesIO()
    scipy.io.mmwrite(b_s, matrix)
    with gzip.open(os.path.join(sample, 'matrix.mtx.gz'), 'wb') as f_out:
        f_out.write(b_s.getvalue())
    write_to_gzip(
        range(len(genes)), os.path.join(sample, 'features.tsv.gz'), lambda i:
        genes[i].decode('UTF-8') + '\t' + gene_names[i].decode('UTF-8') + '\n')
    write_to_gzip(barcodes, os.path.join(sample, 'barcodes.tsv.gz'),
                  lambda i: i.decode('UTF-8') + '\n')


def write_to_gzip(iter, output, line_mod=lambda x: x):
    with gzip.open(output, 'wb') as f_out:
        f_out.writelines([line_mod(i).encode('utf8') for i in iter])


if __name__ == "__main__":
    filename = sys.argv[1]
    data = get_matrix_from_h5(filename)
    writefile(filename, data)
