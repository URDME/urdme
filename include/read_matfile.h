#ifndef READ_MATFILE__H
#define READ_MATFILE__H
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>   /* errno */
extern int errno;
#include <string.h>  /* strerror */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "zlib.h"    /* support for compression (-v7 files) */
#include <time.h>
//**************************************************************************************
//**************************************************************************************
// Need to define some constants
#define mxREAL 10001
#define mxCOMPLEX 10002
//**************************************************************************************
// make sure to set the following environmental variable: SYSTEM_PLATFORM, URDME_VERSION_STRING
//TODO: make this automatic
#define SYSTEM_PLATFORM "MACI64"
#define URDME_VERSION_STRING "1.2b"
//**************************************************************************************
typedef struct{
    short int array_class;
    short int array_flags;
    unsigned int type;
    unsigned int size;
    int is_complex;
    char*name;
    int namesz;
    int ndims;
    int*dims;
    int nvals;
    int nnz;
    size_t* ir;
    int irsz;
    size_t* jc;
    int jcsz;
    double*pr;
    int prsz;
    double*pi;
    int pisz;
} mxArray;
//------------
typedef struct{
    int fd;
    const char*fname;
    short int has_changed;
    short int open_for_writing;
    char header_str[128];
    int num_vars;
    mxArray**vars;
} MATFile;
//**************************************************************************************
int matClose_serialize_array(mxArray*var, unsigned char**ptr);
//----------------------
int matClose_serialize_sparse_array(mxArray*var, unsigned char**ptr);
//----------------------
void matClose(MATFile*f);
//----------------------
void mxDestroyArray(mxArray* var);
//----------------------
MATFile* matOpen(const char*fname,const char* read_write_fg);
//-----------------------------------------------------------------
int matGetNextVariable__read_tag(unsigned int*dtype,unsigned int*dsize,unsigned char*data_ptr);
//-----------------------------------------------------------------
int matGetNextVariable_int32_array(int*outsz,size_t**out,int dsz,unsigned char*data_ptr);
//-----------------------------------------------------------------
int matGetNextVariable_int_array_w_sz(int*outsz,int**out,int dsz,unsigned char*data_ptr,int nvals);
//-----------------------------------------------------------------
int matGetNextVariable_uint_array_w_sz(int*outsz,unsigned int**out,int dsz,unsigned char*data_ptr,int nvals);
//-----------------------------------------------------------------
int matGetNextVariable_int_array(int*outsz,int**out,int dsz,unsigned char*data_ptr);
//-----------------------------------------------------------------
int matGetNextVariable_uint_array(int*outsz,unsigned int**out,int dsz,unsigned char*data_ptr);
//-----------------------------------------------------------------
int matGetNextVariable_double_array(int*outsz,double**out,int dsz,unsigned char*data_ptr);
//-----------------------------------------------------------------
void matGetNextVariable_array(mxArray*mat_var,int dsz,unsigned char*data_ptr);
//----------------------
void matGetNextVariable_sparse_array(mxArray*mat_var,int dsz,unsigned char*data_ptr);
//----------------------
mxArray* matGetNextVariable(MATFile*matfile);
//----------------------
mxArray* matGetVariable(MATFile*matfile,const char* var_name);
//----------------------
size_t* mxGetIr(mxArray* var);
//----------------------
size_t* mxGetJc(mxArray* var);
//----------------------
double *mxGetPr(mxArray* var);
//----------------------
int mxGetNumberOfElements(mxArray* var);
//----------------------
//Get maximum nonzero elements for sparse numeric array
int mxGetNzmax(mxArray* var);
//----------------------
//Get row dimension
int mxGetM(mxArray* var);
//----------------------
//Get column dimension
int mxGetN(mxArray* var);
//----------------------
//Gheck is mxArray is empty
int mxIsEmpty(mxArray* var);
//----------------------
double mxGetScalar(mxArray* var);
//----------------------
// mxCreateDoubleMatrix() 2-D, double-precision, floating-point array initialized to 0
mxArray* mxCreateDoubleMatrix(int rows,int cols,int complexity);
//----------------------
// matPutVariable() Array to MAT-file
void matPutVariable(MATFile*file,const char*name,mxArray*data);
//**************************************************************************************
//**************************************************************************************
#endif
