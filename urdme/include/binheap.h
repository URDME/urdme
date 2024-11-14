/* binheap.h - URDME binary heap */

/* J. Cullhed 2008-06-18. */

#ifndef BINHEAP__H
#define BINHEAP__H

void initialize_heap(double *data,int *INDEX,int *INDEX2,int N);
void percolate_down(int node,double *data,int *INDEX,int *INDEX2,int size);
void percolate_up(int node,double *data,int *INDEX,int *INDEX2,int N);
void update(int node,double *data,int *INDEX,int *INDEX2,int N);
int test_heap_prty(double *data,int *INDEX,int N);
void print_heap(double *data,int N);
void test_min_prty(double *data,int N);

#endif /* BINHEAP__H */
