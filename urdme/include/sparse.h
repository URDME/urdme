/* 
 *  
 * B. Drawert 
 * A. Hellander
 *
 */

#ifndef SPARSE__H
#define SPARSE__H

typedef struct{
	
	int nnz;
	int M;
	int N;
	
	/* Column pointer */
	size_t *jc; 
	/* Row pointer */
	size_t *ir;
	int *s;
	
	
}ccs_matrix;

ccs_matrix *create_ccs_matrix(int M, int N, int nnz);
void free_ccs_matrix(ccs_matrix *m);


#endif