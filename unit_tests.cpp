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

  if(argc<2)
    {
      cerr<<"usage: "<<argv[0]<<" mychunk"<<endl;
    }


  
  
  double *testarrayA,*testarrayB,*testarrayC,*bootA,*bootB;
  int mychunk;
  int totalcasesize=177;
  int totalsnpsize=21000;
  int totalgenesize=1000;
  int bootsize=40;
  int snpstart,snpend;
  int genestart,geneend;

  int snpchunk=1;
  int snpchunksize;
  int genechunk=1;
  int genechunksize;
  int i=0;
  int mysnp;
  int mygene;

  
  const char testfile[]="/home/nwk2/mkl_test/Bootqtl/testmat.h5";
  VSLStreamStatePtr stream;
  VSLSSTaskPtr snptask,genetask;
  MKL_INT x_storage=VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT estimate;
  bool bootCols=true;
  int errcode;
  int *bootrows;
  double *bootmatrix,*tquantiles;
  double mat_initial,mat_elapsed;
  double boot_initial,boot_elapsed;
  int arraysize=64;

  mychunk = 0;
  
  mysnp = mychunk/(genechunk+1);
  mygene = mychunk%(genechunk+1);
  
  snpchunksize  = totalsnpsize/snpchunk;
  genechunksize = totalgenesize/genechunk;

  snpstart = mysnp*snpchunksize;
  snpend = snpstart+snpchunksize;
  genestart = mygene*genechunksize;
  geneend = genestart+genechunksize;
  
 



  testarrayA = (double*) mkl_malloc((size_t)totalcasesize*(size_t)snpchunksize*sizeof(double),arraysize);
  testarrayB = (double*) mkl_malloc((size_t)totalcasesize*(size_t)genechunksize*sizeof(double),arraysize);
  testarrayA=TestGenerate(testarrayA,totalcasesize,snpchunksize,0,snpstart,1);
  testarrayB=TestGenerate(testarrayB,totalcasesize,genechunksize,0,genestart,1);
  testarrayC = (double*) mkl_malloc((size_t)genechunksize*(size_t)snpchunksize*(size_t)bootsize*sizeof(double),arraysize);
  bootA = (double*) mkl_malloc((size_t)totalcasesize*(size_t)snpchunksize*sizeof(double),arraysize);
  bootB = (double*) mkl_malloc((size_t)totalcasesize*(size_t)genechunksize*sizeof(double),arraysize);
  bootrows = (int*) mkl_malloc((size_t)totalcasesize*(size_t)bootsize*sizeof(int),arraysize);

  bootrows=MakeBootRows(bootrows,totalcasesize,bootsize);
  //  PrintMat(bootrows,bootsize,totalcasesize);
  //  PrintMat(testarrayA,rowsize,colsize);

  //  PrintMat(testarrayA,totalcasesize,snpchunksize);
  //  PrintMat(testarrayB,totalcasesize,genechunksize);
  mat_elapsed=0;
  boot_elapsed=0;



  for(i=0; i<bootsize; i++)
    {
   
      cout<<"Iteration: "<<i<<endl;
      bootB=DoBootstrap(bootB,testarrayB,totalcasesize,genechunksize,&bootrows[index(i,0,totalcasesize)]);
      bootA=DoBootstrap(bootA,testarrayA,totalcasesize,snpchunksize ,&bootrows[index(i,0,totalcasesize)]);      
  
      boot_initial=dsecnd();

      boot_elapsed+=dsecnd()-boot_initial;
      mat_initial=dsecnd();


      if(bootCols)
	{
	  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,snpchunksize,genechunksize,totalcasesize,1,bootA,snpchunksize,bootB,genechunksize,0,&testarrayC[snpchunksize*genechunksize*i],genechunksize);
	}
      else
	{
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,snpchunksize,genechunksize,totalcasesize,1,bootA,totalcasesize,bootB,genechunksize,0,&testarrayC[snpchunksize*genechunksize*i],genechunksize);
	}
      mat_elapsed+=dsecnd()-mat_initial;
    }
  mkl_free(bootrows);
  cout<<"bootrows free"<<endl;

  mkl_free(testarrayA);
  mkl_free(testarrayB);
  mkl_free(testarrayC);
  mkl_free(bootA);
  mkl_free(bootB);



  cout<<"Boot time: "<<boot_elapsed<<endl;
  cout<<"Mult time: "<<mat_elapsed<<endl;


  return(0);
}
