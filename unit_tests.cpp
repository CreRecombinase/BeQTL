//Tests for functions in mktest.hpp

#include<iostream>
#include "mktest.hpp"
#include "stdlib.h"
#include "mkl.h"
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;

int main(int argc, char* argv[])
{
  /*
  if(argc<2)
    {
      cerr<<"usage: "<<argv[0]<<" snpchunk"<<endl;
    }
  */

  
  
  double *testarrayA,*testarrayB;
  //  int snpchunk;
  int rowsize=3;
  int colsize=3;
  int bootsize=3;
  //  int mysize=1000;

  const char testfile[]="/home/nwk2/mkl_test/Bootqtl/testmat.h5";
  VSLStreamStatePtr stream;
  VSLSSTaskPtr snptask,genetask;
  MKL_INT x_storage=VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT estimage;
  int errcode;
  int *bootrows;
  double *bootmatrix,*tquantiles;

  int colstart,colend;
  
  //  snpchunk = atoi(argv[1]);
  //  colstart = snpchunk*mysize;
  //  colend=colstart+mysize-1;
  
  hid_t file,dataset,filespace,status;
  hsize_t setdims[3]={rowsize,colsize,bootsize};
  file=-1;
      /*
  if(snpchunk==0)
    {
      file = H5Fcreate(testfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      filespace = H5Screate_simple(3,setdims,NULL);
      dataset = H5Dcreate2(file,"cor",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      cout<<"First file access return: "<<file<<endl;

    } else
    {
      while(file<0)
	{
	  file = H5Fopen(testfile,H5F_ACC_RDWR,H5P_DEFAULT);
	  usleep(100);
	}
      dataset = H5Dopen(file,"cor",H5P_DEFAULT);
    }
      */  



  testarrayA = (double*) mkl_malloc(bootsize*rowsize*colsize*sizeof(double),64);


  //  PrintMat(testarrayA,rowsize,colsize);
  MakeBootRows(bootrows,rowsize,bootsize);


  for(int i=0; i<bootsize; i++)
    {
      cout<<"Iteration: "<<i<<endl;
      TestGenerate(&testarrayA[multindex(0,0,i,rowsize,colsize)],rowsize,colsize,colstart,i);
      //PrintMat(&testarrayA[multindex(0,0,i,rowsize,colsize)],rowsize,colsize);
      //DoBootstrap(testarrayB,testarrayA,rowsize,colsize,&bootrows[index(i,0,rowsize)]);
      //DoScale(testarrayB,rowsize,colsize);
      //      WriteMat(testarrayA,0,rowsize-1,0,colsize-1,i,dataset,filespace);
      
    }
  /*
  errcode = H5Sclose(filespace);
  errcode = H5Dclose(dataset);
  errcode = H5Fclose(file);
  */


  //bootmatrix = (double*)mkl_malloc(rowsize*bootsize*sizeof(double),64);
  PrintMat(testarrayA,bootsize,rowsize*colsize);
  tquantiles = (double*)mkl_malloc(2*colsize*rowsize*sizeof(double),64);
  
  
      //      MatSlice(bootmatrix,testarrayA,i,colsize,rowsize,bootsize);
      //      GetSlice(testfile,0,i,0,rowsize,1,bootsize,bootmatrix);
  GetQuantile(testarrayA,rowsize,colsize,bootsize,tquantiles);
  //  PrintMat(tquantiles,rowsize,2);
  PrintMat(tquantiles,rowsize*colsize,2);
  


  return(0);
}

  
