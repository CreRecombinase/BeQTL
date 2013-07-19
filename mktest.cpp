#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl.h"
#include "mktest.hpp"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <sys/file.h>



using namespace std;

int main(int argc, char* argv[])
{
  double *Snpmat=NULL, *Expmat=NULL, *C=NULL;
  double *BootSNP=NULL, *BootExp=NULL;
  
  //annotation parameters
  int *snparray,*genearray;  
  string *snpname,*genename;
  bool *cistrans;
  double mindist= 100000;

  
  //Matrices and temporary vectors for reporting and storing quantiles
  double *tquantile;
  
  


  //file parameters
  char annofilename[]="/scratch/nwk2/hdf5/snpgeneanno.h5";
  char snpgenefilename[]="/scratch/nwk2/hdf5/snpgenemat_BRCAN.h5";
  char quantilefilename[]="/scratch/nwk2/hdf5/quantilemat.h5";
  char readfilelock[]="/scratch/nwk2/hdf5/readlock.txt";
  char writefilelock[]="/scratch/nwk2/hdf5/writelock.txt";

  int wlockf=-1;
  int result=-1;

  // Platform data/case numbers
  int casetotal = 177;
  int snptotal = 906598;
  int genetotal = 20501;

  //variables for profiling
  double s_initial_boot,s_elapsed_boot;
  double s_initial_mult, s_elapsed_mult;
  double s_initial_quant,s_elapsed_quant;
  double s_initial_write,s_elapsed_write;
  int snpchunk;

  snpchunk = atoi(argv[1]);

  int snpsize = 5000;
  int genesize = 5000;



  int snpstart;
  int casestart=0;
  int genestart=0;
  int bsi = 40;

  snpstart =snpchunk*snpsize;
  int casesize = 177;

  //Parameters for matrix multiplication using dgemm
  int alpha = 1.0/((double)casesize), beta = 0.0;  
  int *bootrows;



  cout<<"Initializing large corarray"<<endl;
  C = (double *)mkl_malloc( snpsize*genesize*bsi*sizeof( double ), 64 );
  cout<<"Filling large corarray"<<endl;
  memset(C,0,snpsize*genesize*bsi*sizeof(double));


  cout<<"Reading in matrices!"<<endl;
  ReadMatrix(Expmat,casetotal,genetotal,casestart,casesize,genestart,genesize,snpgenefilename,"genes",readfilelock);  
  ReadMatrix(Snpmat,casetotal,snptotal,casestart,casesize,snpstart,snpsize,snpgenefilename,"snps",readfilelock);
  cout<<"Matrices read!"<<endl;
  
  //Bootstrapping loop
  //Initialize Bootstrap Matrices


  BootSNP= (double *) mkl_malloc(casesize*snpsize*sizeof(double),64);
  BootExp= (double *) mkl_malloc(casesize*genesize*sizeof(double),64);
  
  MakeBootRows(bootrows,casesize,bsi);

  s_elapsed_boot=0;
  s_elapsed_mult=0;
  //  s_elapsed_write=0;
  for( int q=0; q<bsi; q++){
    cout<<"BSI:"<<q<<endl;
    s_initial_boot=dsecnd();
    DoBootstrap(BootSNP,Snpmat,casesize,snpsize,&bootrows[index(q,0,casesize)]);
    DoBootstrap(BootExp,Expmat,casesize,genesize,&bootrows[index(q,0,casesize)]);
    DoScale(BootSNP,casesize,snpsize);
    DoScale(BootExp,casesize,genesize);
    s_elapsed_boot+=(dsecnd()-s_initial_boot);
    s_initial_mult=dsecnd();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, snpsize, genesize, casesize, (1.0/(double) casesize), BootSNP, snpsize, BootExp, genesize, beta, &C[multindex(0,0,q,snpsize,genesize)], genesize);
    s_elapsed_mult+=(dsecnd()-s_initial_mult);
  }
 
  printf ("\n Deallocating memory \n\n");
  mkl_free(Snpmat);
  mkl_free(Expmat);
  mkl_free(BootSNP);
  mkl_free(BootExp);
 
  //Initializing quantile matrix
  tquantile = (double*)mkl_malloc(2*snpsize*genesize*sizeof(double),64);

  s_initial_quant=dsecnd();
  GetQuantile(C,snpsize,genesize,bsi,tquantile);
  s_elapsed_quant=dsecnd()-s_initial_quant;

  s_initial_write=dsecnd();
  WriteMat(tquantile,snpstart,snpsize,snpstart,snpsize,snptotal,genetotal,quantilefilename,snpchunk,writefilelock);
  s_elapsed_write=dsecnd()-s_initial_write;

  cout<<"Bootstrap time: "<<s_elapsed_boot<<endl;
  cout<<"Mult time: "<<s_elapsed_mult<<endl;
  cout<<"Quantile time: "<<s_elapsed_quant<<endl;
  cout<<"Write time: "<<s_elapsed_write<<endl;
  return (0);
}
