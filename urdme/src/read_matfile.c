/* Brian Drawert */
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>   /* errno */
extern int errno;
#include <string.h>  /* strerror */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
/*#include "zlib.h"  // support for compression (-v7 files) */
#include <time.h>
#include "read_matfile.h"
//**************************************************************************************
//**************************************************************************************
int matClose_serialize_array(mxArray*var, unsigned char**ptr){
    //-----------
    int total_size=0;
    unsigned char*ptr2;
    int c=0;
    //-----------
    //array flags
    unsigned int array_flags_tag[2] = {6,8};
    unsigned char array_flags[8] = {6,0,0,0,0,0,0,0};  //mxDOUBLE_CLASS
    total_size+=16;
    //dimensions array
    unsigned int dimensions_array_tag[2] = {4,8}; // only 2D arrays
    int dimensions_array[2];
    dimensions_array[0] = mxGetM(var);
    dimensions_array[1] = mxGetN(var);
    total_size+=16;
    //array name
    unsigned int array_name_tag[2];
    array_name_tag[0]=1;
    array_name_tag[1]=(unsigned int) var->namesz;
    int array_name_padding=0;
    while( (var->namesz +  array_name_padding)%8!=0){ array_name_padding++; } 
    total_size = total_size+8+var->namesz+array_name_padding;
    //pr
    if(var->prsz==0){ fprintf(stderr,"matClose_serialize_array(): var->prsz = 0\n");exit(0); }
    unsigned int pr_tag[2];
    pr_tag[0]=9; //miDOUBLE
    pr_tag[1]=(unsigned int) var->prsz*8;
    total_size = total_size+8+var->prsz*8;
    //pi
    //TODO: if(var-prsz>0){}
    //---------------------------------
    *ptr=(unsigned char*)calloc(total_size,sizeof(unsigned char));
    ptr2=*ptr;
    //---------------------------------
    memcpy(&ptr2[c],array_flags_tag,8*sizeof(unsigned char));  c+=8;
    //printf("matClose_serialize_array(): %x\t\tarray_flags_tag\n",array_flags_tag);
    memcpy(&ptr2[c],array_flags,8*sizeof(unsigned char));  c+=8;
    //printf("matClose_serialize_array(): %x\t\tarray_flags\n",array_flags);
    memcpy(&ptr2[c],dimensions_array_tag,8*sizeof(unsigned char));  c+=8;
    //printf("matClose_serialize_array(): %x\t\tdimensions_array_tag\n",dimensions_array_tag);
    memcpy(&ptr2[c],dimensions_array,8*sizeof(unsigned char));  c+=8;
    //printf("matClose_serialize_array(): %x\t\tdimensions_array\n",dimensions_array);
    memcpy(&ptr2[c],array_name_tag,8*sizeof(unsigned char));  c+=8;
    //printf("matClose_serialize_array(): %x\t\tarray_name_tag\n",array_name_tag);
    memcpy(&ptr2[c],var->name,var->namesz*sizeof(unsigned char));  c=c + var->namesz + array_name_padding;
    //printf("matClose_serialize_array(): %x\t\tarray_name (%i bytes + %i padding)\n",var->name,var->namesz,array_name_padding);
    memcpy(&ptr2[c],pr_tag,8*sizeof(unsigned char));  c+=8;
    memcpy(&ptr2[c],var->pr,var->prsz*8*sizeof(unsigned char));
    //printf("matClose_serialize_array(): prsz=%i wrote %i bytes\n",var->prsz,var->prsz*8*sizeof(unsigned char));
    //---------------------------------
    return total_size;

}
//----------------------
int matClose_serialize_sparse_array(mxArray*var, unsigned char**ptr){
    //TODO: implement
    //unsigned int array_flags_tag[2] = {6,8};
    //unsigned char array_flags[8] = {0,0,0,5,0,0,0,0};  //mxSPARSE_CLASS
    return 0;
}
//**************************************************************************************
void matClose(MATFile*f){
    int fd,i,j,dsz,wr,wsz,ptr_p;
    unsigned char tag[8];
    int fourteen=14;
    unsigned char*ptr;
    close(f->fd);
    if(f->open_for_writing && f->has_changed){ // write out file
        //open file
        //fd = open(f->fname,O_CREAT|O_WRONLY); //open the file for writing, and remove any previous file
        fd = open(f->fname,O_CREAT | O_WRONLY , 
                S_IRUSR | S_IWUSR); //open the file for writing, and remove any previous file
        //write header
        if(write(fd,f->header_str,128)!=128){ fprintf(stderr,"write error, header\n");exit(0); }
        // loop for f->vars, write out each
        for(i=0;i<f->num_vars;i++){
            //serialize data
            dsz = matClose_serialize_array( f->vars[i], &ptr );
            //fprintf(stderr,"%i = matClose_serialize_array();\n",dsz);
            if(dsz==0){ fprintf(stderr,"matClose(): error writing variable '%s', matClose_serialize_array() return 0 bytes\n",f->vars[i]->name);exit(0); }
            //make the tag
            for(j=0;j<8;j++){ tag[j]=0; } // clean it up
            memcpy(tag,&fourteen,4*sizeof(char));  //type=14 is a matrix type
            memcpy(&tag[4],&dsz,4*sizeof(char));   //size 
            //write out tag
            if(write(fd,&tag,8) != 8){ fprintf(stderr,"write error, tag for var %i\n",i);exit(0); }
            //write out data
            wr=dsz;
            ptr_p=0;
            //fprintf(stderr,"\twriting out data for var %i\n",i);
            do{
                wsz = write(fd,&ptr[ptr_p],wr);
                //fprintf(stderr,"\twrote %i bytes\n",wsz);
                if(wsz==0){ fprintf(stderr,"write error, data for var %i\n",i);exit(0); }
                wr-=wsz;
                ptr_p+=wsz;
            }while(wr>0);
            free(ptr);
        }
        free(f->vars);
        //close file
    }
    free(f);
}
//----------------------
void mxDestroyArray(mxArray* var){
    
    if (var != NULL){
      if(var->namesz>0) free(var->name);
      if(var->ndims>0) free(var->dims);
      if(var->irsz>0) free(var->ir);
      if(var->jcsz>0) free(var->jc);
      if(var->prsz>0) free(var->pr);
      if(var->pisz>0) free(var->pi);
      free(var);
    }
}
//----------------------
MATFile* matOpen(const char*fname,const char* read_write_fg){
    MATFile* ret;
    //char header_str[128];
    long file_len,cur_pos;
    int ndx;
    short int endian,version;
    mxArray*tmp_var=NULL;
    time_t t;
    //-----
    //fprintf(stderr,"matOpen(%s,%s)\n",fname,read_write_fg);
    //-----
    ret=(MATFile*)calloc(1,sizeof(MATFile));
    ret->fname = fname;
    if(strcmp(read_write_fg,"r")==0){
        ret->fd = open(fname,O_RDONLY);
        if(ret->fd<0){
            fprintf(stderr,"Can not open file '%s' : %s\n",fname,strerror(errno));
            return NULL;
        }
        // check that we are using v5.0 MAT files ( -v6 & -v7 )
        if( read(ret->fd,&ret->header_str,128*sizeof(char)) != 128){ fprintf(stderr,"ERROR: short read on header\n");exit(0); }
        if( strncmp( ret->header_str , "MATLAB 5.0 MAT-file" , 19 )!=0){ fprintf(stderr,"ERROR: can not read this version of MAT file, use '-v6' or '-v7' with save() command\n\n"); exit(0); }
        return ret;
        //-----------------
    }else if(strcmp(read_write_fg,"w")==0){
        ret->open_for_writing=1;
        // check if the file exists
        ret->fd = open(fname,O_RDONLY);
        if(ret->fd<0){ // it does not
            memset(ret->header_str,' ',128);
            t = time(NULL);
            ndx=snprintf(ret->header_str,116,"MATLAB 5.0 MAT-file, Platform: %s, Created By: URDME-libmat v%s on %s", SYSTEM_PLATFORM, URDME_VERSION_STRING,ctime(&t));
            //ret->header_str[115] = '\0';
            ret->header_str[ndx] = ' ';
            endian = 0x4d49;
            version = 0x0100;
            memcpy(&ret->header_str[124],&version,2);
            memcpy(&ret->header_str[126],&endian,2);
            //ret->header_str, "MATLAB 5.0 MAT-file, Platform: URDME_LIBMAT Created on: TODO                                                                  IM");
            //if(strlen(ret->header_str) != 128){ fprintf(stderr,"header string is not strlen() of 128!, strlen=%lu\n",strlen(ret->header_str));exit(0); }
            return ret;
        }else{ // it does
            if( read(ret->fd,&ret->header_str,128*sizeof(char)) != 128){ fprintf(stderr,"ERROR: short read on header\n");exit(0); }
            // read all the variables and store them in the MATFile struct
            file_len  = (long) lseek(ret->fd,0,SEEK_END);
            lseek(ret->fd, 128, SEEK_SET);
            do{
                cur_pos =  (long) lseek(ret->fd,0,SEEK_CUR);
                if(cur_pos>=file_len){ break; }
                tmp_var = matGetNextVariable(ret);
                matPutVariable(ret, tmp_var->name, tmp_var );
            }while(tmp_var!=NULL);
            ret->has_changed=0;
            return ret;
        }
    }else{
        fprintf(stderr,"matOpen: flag is not 'r' or 'w'\n");exit(0);
    }
}
//***********************************************************************************************

int matGetNextVariable__read_tag(unsigned int*dtype,unsigned int*dsize,unsigned char*data_ptr){
    int tag[2];
    memcpy(tag,data_ptr,8);
    *dtype = tag[0] & 0x000000ff;
    if ( tag[0] & 0xffff0000 ) { /* Data is packed in the tag */
        *dsize = (tag[0] >> 16);
        return 4;
    }else{
        *dsize=tag[1];
        return 8;
    }
}
//-----------------------------------------------------------------
// use to get IR and JC arrays, note the output of these arrays is size_t, not int32
int matGetNextVariable_int32_array(int*outsz,size_t**out,int dsz,unsigned char*data_ptr){
    unsigned int dtype,dsize,nvals;
    int utmp;
    int i,offset;
    size_t*out2;
    offset = matGetNextVariable__read_tag(&dtype,&dsize,data_ptr);
    if(dtype!=5){fprintf(stderr,"Error: array type not miINT32, type %i\n",dtype);exit(0);}
    nvals = dsize/4;
    //printf("print_int32_array(): %i values\n",nvals);
    *out = (size_t*)malloc(nvals*sizeof(size_t));
    out2=*out;
    if(nvals>dsz){fprintf(stderr,"Error: array size mis-match. dsz=%i\n",dsz);exit(0);}
    for(i=0;i<nvals;i++){
        utmp=0;
        memcpy(&utmp,&data_ptr[offset+4*i],4*sizeof(char));
        out2[i]=(size_t)utmp;    
        //printf("%i ",utmp);
    }
    //printf("\n");
    *outsz=nvals;
    return offset+dsize;
}
//-----------------------------------------------------------------
//  double,uint16,uint8,int8  :: return double array
//  uint32  :: return int array
//  probably need to have it go both ways
//-----------------------------------------------------------------
// take next data segment as a int "out", returns array size
int matGetNextVariable_int_array_w_sz(int*outsz,int**out,int dsz,unsigned char*data_ptr,int nvals){
    unsigned int dtype,dsize,tsize;
    int itmp;
    int*out2;
    double dtmp;
    int i,offset;
    offset = matGetNextVariable__read_tag(&dtype,&dsize,data_ptr);
    //printf("matGetNextVariable_int_array(): dtype=%i  addy=%d\n",dtype,data_ptr);
    if(dtype==9){ // type convert from double
        tsize = dsize/nvals;
        *out = (int*)malloc(nvals*sizeof(int));
        out2=*out;
        for(i=0;i<nvals;i++){
            memcpy(&dtmp,&data_ptr[offset+i*tsize],tsize*sizeof(char));
            out2[i] =(int)dtmp;
        }
        *outsz=nvals;
        return offset+dsize;
    }else{ //direct copy
        tsize=dsize/nvals;
        printf("dtype=%i tsize=%i nvals=%i\n",dtype,tsize,nvals);
        *out = (int*)malloc(nvals*sizeof(int));
        out2=*out;
        for(i=0;i<nvals;i++){
            itmp=0;
            memcpy(&itmp,&data_ptr[offset+i*tsize],tsize*sizeof(char));
            out2[i]=itmp;
        }
        
        *outsz=nvals;
        return offset+dsize;
    }
}
//-----------------------------------------------------------------
// take next data segment as a uint "out", returns array size
int matGetNextVariable_uint_array_w_sz(int*outsz,unsigned int**out,int dsz,unsigned char*data_ptr,int nvals){
    // I think we can just cast the uint to an int
    return matGetNextVariable_int_array_w_sz(outsz,(int**)out,dsz,data_ptr,nvals);
}
//-----------------------------------------------------------------
//-----------------------------------------------------------------
//  double,uint16,uint8,int8  :: return double array
//  uint32  :: return int array
//  probably need to have it go both ways
//-----------------------------------------------------------------
// take next data segment as a int "out", returns array size
int matGetNextVariable_int_array(int*outsz,int**out,int dsz,unsigned char*data_ptr){
    unsigned int dtype,dsize,nvals,tsize;
    int itmp;
    int*out2;
    double dtmp;
    int i,offset;
    offset = matGetNextVariable__read_tag(&dtype,&dsize,data_ptr);
    //printf("matGetNextVariable_int_array(): dtype=%i  addy=%d\n",dtype,data_ptr);
    if(dtype==9){ // type convert from double
        tsize=4;
        nvals = dsize/tsize;
        *out = (int*)malloc(nvals*sizeof(int));
        out2=*out;
        for(i=0;i<nvals;i++){
            memcpy(&dtmp,&data_ptr[offset+i*tsize],tsize*sizeof(char));
            out2[i] =(int)dtmp;
        }
        *outsz=nvals;
        return offset+dsize;
    }else if(dtype==5 || dtype==6){ //direct copy
        tsize=2;
        nvals = dsize/tsize;
        printf("dtype=%i tsize=%i nvals=%i\n",dtype,tsize,nvals);
        *out = (int*)malloc(nvals*sizeof(int));
        out2=*out;
        for(i=0;i<nvals;i++){
            itmp=0;
            memcpy(&itmp,&data_ptr[offset+i*tsize],4*sizeof(char));
            out2[i]=(int)itmp;
        }
        
        *outsz=nvals;
        return offset+dsize;
    }else if(dtype==3 || dtype==4){
        tsize=2;
    }else if(dtype==1 || dtype==2){
        tsize=1;
    }else{
        fprintf(stderr,"matGetNextVariable_int_array(), Error: array type not supported, type %i\n",dtype);exit(0);
    }
    nvals = dsize/tsize;
    *out = (int*)malloc(nvals*sizeof(int));
    out2=*out;
    for(i=0;i<nvals;i++){
        out2[i]=(int)data_ptr[offset+i*tsize];
    }
    *outsz=nvals;
    return offset+dsize;
}
//-----------------------------------------------------------------
// take next data segment as a uint "out", returns array size
int matGetNextVariable_uint_array(int*outsz,unsigned int**out,int dsz,unsigned char*data_ptr){
    // I think we can just cast the uint to an int
    return matGetNextVariable_int_array(outsz,(int**)out,dsz,data_ptr);
}
//-----------------------------------------------------------------
// take next data segment as a double "out", returns array size
int matGetNextVariable_double_array(int*outsz,double**out,int dsz,unsigned char*data_ptr){
    unsigned int dtype,dsize,nvals,utmp,tsize;
    //double dtmp;
    int itmp,i,offset;
    double *go_out;
    offset = matGetNextVariable__read_tag(&dtype,&dsize,data_ptr);
    //printf("dtype=%i\n",dtype);
    if(dtype==9){ // direct copy
        nvals = dsize/8;
        *out = (double*)malloc(nvals*sizeof(double));
        go_out=*out;
        memcpy(go_out,&data_ptr[offset],nvals*8*sizeof(char));
        *outsz=nvals;
        return offset+dsize;
    }else if(dtype>=1&&dtype<=6){ // type convert (u)int->double
        if(dtype==5 || dtype==6){ tsize=4;
        }else if(dtype==3 || dtype==4){ tsize=2;
        }else if(dtype==1 || dtype==2){ tsize=1; }else{ fprintf(stderr,"should not get here\n");exit(0); }
        nvals=dsize/tsize;
        *out = (double*)malloc(nvals*sizeof(double));
        go_out=*out;
        if(dtype==1||dtype==3||dtype==5){
            for(i=0;i<nvals;i++){
                itmp=0;
                memcpy(&itmp,&data_ptr[offset+i*tsize],tsize*sizeof(char));
                go_out[i]=(double)itmp;
            }
        }else if(dtype==2||dtype==4||dtype==6){
            for(i=0;i<nvals;i++){
                utmp=0;
                memcpy(&utmp,&data_ptr[offset+i*tsize],tsize*sizeof(char));
                go_out[i]=(double)utmp;
            }
        }
        *outsz=nvals;
        return offset+dsize;
    }else{
        fprintf(stderr,"matGetNextVariable_double_array(), Error: array type not supported, type %i\n",dtype);exit(0);
    }
}
//-----------------------------------------------------------------
void matGetNextVariable_array(mxArray*mat_var,int dsz,unsigned char*data_ptr){
    int offset=0;
    //printf("matGetNextVariable_array()\n");
    //Real Part data //////////////////////////////////////
    offset=matGetNextVariable_double_array(&(mat_var->prsz),&(mat_var->pr),dsz,data_ptr);

    //Complex Part data ////////////////////////////////
    if(mat_var->is_complex){
        while(offset%8 != 0){ offset+=1; }
        if(offset>=dsz){printf("Error: array is complex, but no data left\n");exit(0);}
        matGetNextVariable_double_array(&(mat_var->pisz),&(mat_var->pi),dsz,&data_ptr[offset]);
    }
}
//-----------------------------------------------------------------
void matGetNextVariable_sparse_array(mxArray*mat_var,int dsz,unsigned char*data_ptr){
    int offset=0;
    // ir
    offset=matGetNextVariable_int32_array(&(mat_var->irsz),&(mat_var->ir),dsz,&data_ptr[offset]);
    while(offset%8 != 0){ offset+=1; }
    // jc
    offset+=matGetNextVariable_int32_array(&(mat_var->jcsz),&(mat_var->jc),dsz,&data_ptr[offset]);
    while(offset%8 != 0){ offset+=1; }
    // pr
    //Real Part data //////////////////////////////////////
    offset=matGetNextVariable_double_array(&(mat_var->prsz),&(mat_var->pr),dsz,&data_ptr[offset]);

    //Complex Part data ////////////////////////////////
    if(mat_var->is_complex){
        while(offset%8 != 0){ offset+=1; }
        if(offset>=dsz){printf("Error: array is complex, but no data left\n");exit(0);}
        matGetNextVariable_double_array(&(mat_var->pisz),&(mat_var->pi),dsz,&data_ptr[offset]);
    }
}
//-----------------------------------------------------------------
//-----------------------------------------------------------------
//-----------------------------------------------------------------
mxArray* matGetNextVariable(MATFile*matfile){
    mxArray*mat_var;
    unsigned int dtype,dsize;
    unsigned char *data_ptr;
    unsigned int dim_type,dim_size,stype,ssize,name_type,name_size;
    int dim_tmp,dim_nnz,cnt,i;
    int offset=0;
    //------------------------------
    mat_var=(mxArray*)calloc(1,sizeof(mxArray));
    //------------------------------
    if( ( read(matfile->fd,&dtype,sizeof(unsigned int)) != sizeof(unsigned int) ) ||
    ( read(matfile->fd,&dsize,sizeof(unsigned int)) != sizeof(unsigned int) ) )
    { return NULL; }
    data_ptr=(unsigned char *)malloc(dsize*sizeof(unsigned char));
    if( read(matfile->fd,data_ptr,dsize*sizeof(unsigned char)) != dsize*sizeof(unsigned char) ){
        free(data_ptr);
        return NULL;
    }
    //------------------------------
    mat_var->type=dtype;
    mat_var->size=dsize;
    //------------------------------
    if(dtype==15){fprintf(stderr,"ERROR: don't support compression yet\n");exit(0);}
    if(dtype!=14){fprintf(stderr,"ERROR: don't support non-Matrix types\n");exit(0);}
    //------------------------------//Array Flags
    memcpy(&stype,&data_ptr[0],4*sizeof(char));
    memcpy(&ssize,&data_ptr[4],4*sizeof(char));
    if(stype!=6){fprintf(stderr,"Error: First tag of matrix not miUINT32, type %i\n",stype);exit(0);}
    memcpy(&mat_var->array_class,&data_ptr[8],1*sizeof(char));
    if(mat_var->array_class==5){
        dim_nnz=0;
        memcpy(&dim_nnz,&data_ptr[12],4*sizeof(char));
        mat_var->nnz=dim_nnz;
    }
    memcpy(&mat_var->array_flags,&data_ptr[9],1*sizeof(char));
    if(mat_var->array_flags!=0 && (mat_var->array_flags & 0xc)){
        mat_var->is_complex=1;
    }
    offset = 16;
    //printf("\tarray_flags=%i\n",&mat_var->array_flags);
    //------------------------------ //Dimensions Array
    memcpy(&dim_type,&data_ptr[offset+0],4*sizeof(char));
    memcpy(&dim_size,&data_ptr[offset+4],4*sizeof(char));
    mat_var->ndims = dim_size/4;
    mat_var->dims = (int *)malloc(mat_var->ndims*sizeof(int));
    cnt=0;
    //printf("\tdims: ");
    for(i=0;i<dim_size;i+=4){
        memcpy(&dim_tmp,&data_ptr[offset+8+i],4*sizeof(char));
        mat_var->dims[cnt]=dim_tmp;
        //printf(" %i",mat_var->dims[cnt]);
        cnt++;
    }
    //printf("\n");
    offset += 8 + dim_size;
    while(offset%8 != 0){ offset+=1; }
    //------------------------------ //Array Name
    offset += matGetNextVariable__read_tag(&name_type,&name_size,&data_ptr[offset]);
    mat_var->name = (char *)malloc((name_size+1)*sizeof(char));
    mat_var->namesz=name_size;
    memcpy(mat_var->name,&data_ptr[offset],name_size*sizeof(char));
    mat_var->name[name_size]='\0';
    offset += name_size;
    while(offset%8 != 0){ offset+=1; }
    //printf("\tName: %s\n",mat_var->name);
    //------------------------------
    if(mat_var->array_class==6){       matGetNextVariable_array(mat_var,dsize-offset,&data_ptr[offset]);  } //returns values of type double
    else if(mat_var->array_class==5){  matGetNextVariable_sparse_array(mat_var,dsize-offset,&data_ptr[offset]);  } // not working!
    //else if(array_class==2){  print_structure(dsize-offset,&data_ptr[offset],ndims,dim_arr);  }
    //else if(array_class==1){  print_cell_array(dsize-offset,&data_ptr[offset],ndims,dim_arr);  }
    else{
        //printf(mat_var->name);
        /* AH: THIS CHECK MAKES IT EXIT IF THERE ARE STRINGS IN THE FILE, EVEN IF WE ARE NOT TRYING TO READ THEM. 
           I HAVE COMMENTED OUT THE CHECK, BUT NOW THERE IS A RISK OF SUPPLYING AN UNSOPPRTED STRUCTURE AND TRYING
           TO READ IT. */
         //fprintf(stderr,"ERROR: can not handle data type (array_class=%i)\n",mat_var->array_class);exit(0);
    }
    //------------------------------
    free(data_ptr);
    return mat_var;
}
//**************************************************************************************
mxArray* matGetVariable(MATFile*matfile,const char* var_name){
    mxArray*mat_var;
    long file_len,cur_pos;
    // Start at beginning
    file_len  = (long) lseek(matfile->fd,0,SEEK_END);
    lseek(matfile->fd, 128, SEEK_SET);
    do{
        cur_pos =  (long) lseek(matfile->fd,0,SEEK_CUR);
        if(cur_pos >= file_len){
            //printf("matGetNextVariable() got EOF\n");
            return NULL;
        }
        //printf("matGetNextVariable()  %ld / %ld \n", cur_pos , file_len );
        mat_var = matGetNextVariable(matfile);
        //printf("\tmatGetNextVariable: found %s\n",mat_var->name);
        if(strcmp(mat_var->name,var_name)==0 ){
            return mat_var;
        }
        
        //printf("\tmxDestroyArray(%s)\n",mat_var->name);
        mxDestroyArray(mat_var);
    }while(mat_var!=NULL);
    //fprintf(stderr,"matGetVariable(): error\n");
    //exit(0);
    return NULL;
}
//----------------------
size_t* mxGetIr(mxArray* var){
    return var->ir;
}
//----------------------
//unsigned int* mxGetIr(mxArray* var){
//    return var->ir;
//}
//----------------------
size_t* mxGetJc(mxArray* var){
    return var->jc;
}
//----------------------
//unsigned int* mxGetJc(mxArray* var){
//    return var->jc;
//}
//----------------------
double *mxGetPr(mxArray* var){
    return var->pr;
}
//----------------------
int mxGetNumberOfElements(mxArray* var){
    return var->prsz;
}
//----------------------
//Get maximum nonzero elements for sparse numeric array
int mxGetNzmax(mxArray* var){
    return var->nnz;
}
//----------------------
//Get row dimension
int mxGetM(mxArray* var){
    return var->dims[0];
}
//----------------------
//Get column dimension
int mxGetN(mxArray* var){
    return var->dims[1];
}
//----------------------
double mxGetScalar(mxArray* var){
    if(var->prsz>0){
        return var->pr[0];
    }
    return -1;
}

int mxIsEmpty(mxArray *var){
    
    if (var == NULL) {
        return 1;
    }
    
    int i;
    if (var->ndims == 0) {
        return 1;
    }
    for (i = 0; i < var->ndims; i++) {
        if (var->dims[i]==0) {
            return 1;
        }
    }
    return 0;
}

//----------------------
// mxCreateDoubleMatrix() 2-D, double-precision, floating-point array initialized to 0
mxArray* mxCreateDoubleMatrix(int rows,int cols,int complexity){
    mxArray*mat_var;
    //---------------
    if(complexity!=mxREAL){
        fprintf(stderr,"mxCreateDoubleMatrix(); mxCOMPLEX type not implimented yet\n");
        exit(0);
    }
    //---------------
    mat_var=(mxArray*)calloc(1,sizeof(mxArray));
    //---------------
    mat_var->type = 14;  //double matrix
    mat_var->size = -1; //this gets set in matPutVariable()
    mat_var->nnz = rows*cols;
    mat_var->ndims = 2;
    mat_var->dims = (int *)malloc(mat_var->ndims*sizeof(int));
    mat_var->dims[0] = rows;
    mat_var->dims[1] = cols;
    //mat_var->namesz = ;  //this gets set in matPutVariable()
    //mat_var->name = (char *)malloc(name_size*sizeof(char));
    mat_var->array_class=6;
    //---------------
    mat_var->pr = (double*)calloc(rows*cols,sizeof(double));
    mat_var->prsz = rows*cols;
    //---------------
    return mat_var;
}
//----------------------
// matPutVariable() Array to MAT-file
void matPutVariable(MATFile*file,const char*name,mxArray*data){
    int i;
    mxArray*tmp_var;
    mxArray**new_vars;
    mxArray**old_vars;
    //fprintf(stderr,"matPutVariable not implimented yet\n");
    //exit(0);
    data->namesz=strlen(name);
    data->name = (char*)malloc(data->namesz*sizeof(char));
    strcpy(data->name,name);
    // iterate through all the variables in the MATFile, if any of the names match, replace that entry
    //fprintf(stderr,"\there1\n");
    for(i=0;i<file->num_vars;i++){
        tmp_var = file->vars[i];
        if( strcmp(tmp_var->name,data->name)==0 ){
            file->vars[i] = data;
            mxDestroyArray(tmp_var);
            return;
        }
    }
    // otherwise allocate a new array of size file->num_vars+1, copy all existing elements over, add 'data' to end
    //fprintf(stderr,"\tmalloc new_vars\n");
    new_vars = (mxArray**)malloc((file->num_vars+1)*sizeof(mxArray*));
    old_vars = file->vars;
    for(i=0;i<file->num_vars;i++){
        new_vars[i] = old_vars[i];
    }
    new_vars[i] = data;
    file->vars = new_vars;
    file->num_vars = file->num_vars+1;
    file->has_changed=1;
    free(old_vars);
}
//**************************************************************************************
//**************************************************************************************
