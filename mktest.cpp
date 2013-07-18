#include <iostream>
#include <stdio.h>
#include "mkl.h"
#include "mktest.hpp"
#include "hdf5.h"
#include "hdf5_hl.h"



using namespace std;

int main()
{
  double *Snpmat=NULL, *Expmat=NULL, *C=NULL;
  double *BootSNP=NULL, *BootExp=NULL;
  
  //annotation parameters
  int *snparray,*genearray;  
  string *snpname,*genename;
  bool *cistrans;
  double mindist= 100000;

  
  //Matrices and temporary vectors for reporting and storing quantiles
  double *lowquantiles,*highquantiles;
  double *tquantile,*bootmatrix;
  
  


  //file parameters
  char annofilename[]="/scratch/nwk2/hdf5/snpgeneanno.h5";
  char corfilename[]="/scratch/nwk2/hdf5/testcor.h5";
  char snpgenefilename[]="/scratch/nwk2/hdf5/snpgenemat_BRCAN.h5";
  char quantilefilename[]="/scratch/nwk2/hdf5/quantilemat.h5";

  // Platform data/case numbers
  int casetotal = 177;
  int snptotal = 906598;
  int genetotal = 20501;

  //variables for profiling
  double s_initial_boot,s_elapsed_boot;
  double s_initial_mult, s_elapsed_mult;
  double s_initial_write,s_elapsed_write;

  int mysnpstart=0;
  int mysnpend=100;
  int mycasestart=0;
  int mycaseend=176;
  int mygenestart=0;
  int mygeneend=2050;
  int bsi = 3;
  
  int snpsize = (mysnpend-mysnpstart)+1;
  int genesize= (mygeneend-mygenestart)+1;
  int casesize = (mycaseend-mycasestart)+1;

  //Parameters for matrix multiplication using dgemm
  int alpha = 1.0/((double)casesize), beta = 0.0;  

  hid_t file, dataset,filespace,status;
  hsize_t cordims[3]={snptotal,genetotal,bsi};

  int *bootrows;

  

  C = (double *)mkl_malloc( snpsize*genesize*sizeof( double ), 64 );
  for (int i = 0; i < (snpsize*genesize); i++) {
    C[i] = 0.0;
  }
  ReadMatrix(Expmat,casetotal,genetotal,mycasestart,mycaseend,mygenestart,mygeneend,snpgenefilename,"genes");  
  ReadMatrix(Snpmat,casetotal,snptotal,mycasestart,mycaseend,mysnpstart,mysnpend,snpgenefilename,"snps");
  //  PrintMat(Expmat,casesize,genesize);
  //  PrintMat(Snpmat,casesize,snpsize);
  

  
  file = H5Fcreate(corfilename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  filespace = H5Screate_simple(3,cordims,NULL);

  dataset = H5Dcreate2(file,"cor",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  

  //Bootstrapping loop
  //Initialize Bootstrap Matrices


  BootSNP= (double *) mkl_malloc(casesize*snpsize*sizeof(double),64);
  BootExp= (double *) mkl_malloc(casesize*genesize*sizeof(double),64);
  
  MakeBootRows(bootrows,casesize,bsi);



  s_elapsed_boot=0;
  s_elapsed_mult=0;
  s_elapsed_write=0;
  for( int q=0; q<bsi; q++){
    cout<<"BSI:"<<q<<endl;
    s_initial_boot=dsecnd();
    DoBootstrap(BootSNP,Snpmat,casesize,snpsize,&bootrows[index(q,0,casesize)]);
    DoBootstrap(BootExp,Expmat,casesize,genesize,&bootrows[index(q,0,casesize)]);
    
    DoScale(BootSNP,casesize,snpsize);
    DoScale(BootExp,casesize,genesize);
    s_elapsed_boot+=(dsecnd()-s_initial_boot);
    s_initial_mult=dsecnd();
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, snpsize, genesize, casesize, (1.0/(double) casesize), BootSNP, snpsize, BootExp, genesize, beta, C, genesize);
    s_elapsed_mult+=(dsecnd()-s_initial_mult);
    //printf ("\n Computations completed.\n\n");
    
    s_initial_write=dsecnd();
    WriteMat(C,mysnpstart,mysnpend,mygenestart,mygeneend,q,dataset,filespace);
    s_elapsed_write+=(dsecnd()-s_initial_write);
  }

  status = H5Dclose(dataset);
  status = H5Sclose(filespace);
  status = H5Fclose(file);
  cout<<"Bootstrap time: "<<s_elapsed_boot<<endl;
  cout<<"Mult time: "<<s_elapsed_mult<<endl;
  cout<<"Write time: "<<s_elapsed_write<<endl;
  
  
  //printf ("\n Deallocating memory \n\n");
  mkl_free(Snpmat);
  mkl_free(Expmat);
  mkl_free(BootSNP);
  mkl_free(BootExp);
  mkl_free(C);

  
  cordims[2]=2;

  file = H5Fcreate(quantilefilename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  filespace = H5Screate_simple(3,cordims,NULL);
  dataset = H5Dcreate2(file,"quantiles",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  lowquantiles = (double*)mkl_malloc(snpsize*genesize*sizeof(double),64);
  cout<<"Low quantiles made"<<endl;
  highquantiles = (double*)mkl_malloc(snpsize*genesize*sizeof(double),64);
  tquantile = (double*)mkl_malloc(snpsize*sizeof(double),64);


  bootmatrix = (double *) mkl_malloc(snpsize*bsi*sizeof(double),64);
  for(int i=0; i<genesize; i++)
    {
      GetSlice(corfilename,mysnpstart,i,0,snpsize,1,bsi,bootmatrix);
      GetQuantile(bootmatrix,snpsize,bsi,tquantile);
      cblas_dcopy(snpsize,tquantile,2,&lowquantiles[index(0,i,genesize)],genesize);
      cblas_dcopy(snpsize,&tquantile[1],2,&highquantiles[index(1,i,genesize)],genesize);
    }


  cout<<"Writing lowquantile"<<endl;
  WriteMat(lowquantiles,mysnpstart,mysnpend,mygenestart,mygeneend,0,dataset,filespace);
  cout<<"Writing highquantile"<<endl;
  WriteMat(highquantiles,mysnpstart,mysnpend,mygenestart,mygeneend,1,dataset,filespace);
  
  status = H5Dclose(dataset);
  status = H5Sclose(filespace);
  status = H5Fclose(file);  

  
  return 0;
}
