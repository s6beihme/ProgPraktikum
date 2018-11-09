#include <stdio.h>
#include <stdlib.h>

//a)

struct _csr_matrix {
	int row_count;
	int col_count;
	int nnz;
	double* entr;
	int* col_ind;
	int* new_row_ind;
};

typedef struct _csr_matrix* CsrMatrix;

//function to create matrix in csr format
//arguments: 1. number of rows in matrix, 2. number of columns in matrix,
//3. number of non-zero entries in matrix, 4. matrix parameter to be created
int csr_create(int rows, int columns, int nnz, CsrMatrix* mat) {
	if (rows < 0 || columns < 0 || nnz < 0 || mat == NULL) {
		printf("csr_create failed due to invalid input");
		return -1;
	}
	(*mat) = malloc(sizeof(CsrMatrix));
	if ((*mat) == NULL) return -1;

	(*mat)->row_count = rows;
	(*mat)->col_count = columns;
	(*mat)->nnz = nnz;
	(*mat)->entr = malloc(nnz * sizeof(double));
	if ((*mat)->entr == NULL) {
		free((*mat));
		return -1;
	}
	(*mat)->col_ind = malloc(nnz * sizeof(int));
	if ((*mat)->col_ind == NULL) {
		free((*mat)->entr);
		free((*mat));
		return -1;
	}
	(*mat)->new_row_ind = malloc(rows * sizeof(int));
	if ((*mat)->new_row_ind == NULL) {
		free((*mat)->entr);
		free((*mat)->col_ind);
		free((*mat));
		return -1;
	}
	return 0;
}

//function to free space taken up by matrix of type CsrMatrix
void csr_free(CsrMatrix* mat) {
	if (mat == NULL || (*mat) == NULL) return;

	free((*mat)->entr);
	free((*mat)->col_ind);
	free((*mat)->new_row_ind);
	free((*mat));
}

//function to initialize matrix of CsrMAtrix type
//parameters:
	//mat: matrix to initialize, 
	//vals: array of length nnz containing values of new entries
	//row_inds: array of length nnz containing row indices of new entries. arranged in ascending order
	//col_inds: array of length nnz containing column indices of new entries. arranged in ascending order per row
	//nnz: number of values to be set
void csr_assemble(CsrMatrix mat, double* vals, int* row_inds, int* col_inds, int nnz) {
	if (mat == NULL || vals == NULL || row_inds == NULL || col_inds == NULL || nnz < 0 || nnz!=mat->nnz) {
		printf("csr_assemble failed due to invalid input");
		return;
	}

	//copy vals and col_inds into matrix, so you can handle parameter arrays and matrix seperately after function
	for (int i = 0; i < nnz; i++) {
		mat->entr[0] = vals[0];
		mat->col_ind[0] = col_inds[0];
	}

	mat->new_row_ind[0] = 0;

	//variable to keep track of where in new_row_inds you currently are
	int cur_pos = 0;

	//fill up new_row_inds untill you reach first no-zero entry
	for (int i = 0; i < row_inds[0]; i++) {
		cur_pos++;
		mat->new_row_ind[cur_pos] = 0;
	}
	

	//create indices of new rows for rest of entries
	for (int i = 1; i < nnz; i++) {
		if (row_inds[i] > row_inds[i - 1]) {
			for (int j = 0; j < (row_inds[i] - row_inds[i - 1]); j++) {
				cur_pos++;
				mat->new_row_ind[cur_pos] = i;
			}
		}
	}
	cur_pos++;
	mat->new_row_ind[cur_pos] = nnz;
}

int main() {
	double v[19] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	int r[19] =    { 0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,4,4,5,5};
	int c[19] =    { 0,1,3,0,1,3,1,2,3,4,3,5,1,2,3,4,5,2,5};
	CsrMatrix m;
	csr_create(6, 6, 19, &m);
	csr_assemble(m, v, r, c, 19);
	for (int i = 0; i < 7; i++) {
		printf("%i, ",m->new_row_ind[i]);
	}

	return 0;
}