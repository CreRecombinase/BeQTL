//Tests for functions in mktest.hpp

#include<iostream>
#include "mktest.hpp"
#include "mkl.h"
//#include "hdf5.h"
//#include "hdf5_hl.h"

int main()
{
  double *testarrayA,*testarrayB;
  int rowsize=4;
  int colsize=5;
  int bootsize=4;
  const char testfile[]="/home/nwk2/mkl_test/Bootqtl/testmat.h5";
  VSLStreamStatePtr stream;
  VSLSSTaskPtr snptask,genetask;
  MKL_INT x_storage=VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT estimage;
  int errcode;
  int *bootrows;
  double *bootmatrix,*tquantiles;
  
  hid_t file,dataset,filespace,status;
  hsize_t setdims[3]={rowsize,colsize,bootsize};
  file = H5Fcreate(testfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  filespace = H5Screate_simple(3,setdims,NULL);
  dataset = H5Dcreate2(file,"cor",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);


  testarrayB = (double*) mkl_malloc(rowsize*colsize*sizeof(double),64);

  TestGenerate(testarrayA,rowsize,colsize,0);
  PrintMat(testarrayA,rowsize,colsize);
  MakeBootRows(bootrows,rowsize,bootsize);


  for(int i=0; i<bootsize; i++)
    {
      //DoBootstrap(testarrayB,testarrayA,rowsize,colsize,&bootrows[index(i,0,rowsize)]);
      //DoScale(testarrayB,rowsize,colsize);
      WriteMat(testarrayA,0,rowsize-1,0,colsize-1,i,dataset,filespace);
      
    }
  bootmatrix = (double*)mkl_malloc(rowsize*bootsize*sizeof(double),64);
  tquantiles = (double*)mkl_malloc(rowsize*sizeof(double),64);
  
  
  for(int i=0; i< colsize; i++)
    {
      GetSlice(testfile,0,i,0,rowsize,1,bootsize,bootmatrix);
      GetQuantile(bootmatrix,rowsize,bootsize,tquantiles);
      PrintMat(tquantiles,rowsize,2);
      PrintMat(bootmatrix,rowsize,colsize);
    }

  return(0);
}

  
