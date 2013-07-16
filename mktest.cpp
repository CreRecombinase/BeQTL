/* C source code is found in dgemm_example.c */


#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "mktest.hpp"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <time.h>





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
  double *quantiles;
  double *tquantile,*bootmatrix;
  
  


  //file parameters
  char annofilename[]="/scratch/nwk2/hdf5/snpgeneanno.h5";
  char corfilename[]="/tmp/nwk2/testcor.h5";
  char snpgenefilename[]="/scratch/nwk2/hdf5/snpgenemat_BRCAN.h5";
  
  // Platform data/case numbers
  int casetotal = 177;
  int snptotal = 906598;
  int genetotal = 20501;

  //variables for profiling
  double s_initial_boot,s_elapsed_boot;
  double s_initial_mult, s_elapsed_mult;
  double s_initial_write,s_elapsed_write;

  int mysnpstart=0;
  int mysnpend=10;
  int mycasestart=0;
  int mycaseend=5;
  int mygenestart=0;
  int mygeneend=5;
  int bsi = 5;
  
  int snpsize = (mysnpend-mysnpstart)+1;
  int genesize= (mygeneend-mygenestart)+1;
  int casesize = (mycaseend-mycasestart)+1;

  //Parameters for matrix multiplication using dgemm
  int alpha = 1.0/((double)casesize), beta = 0.0;  

  VSLStreamStatePtr stream;
  VSLSSTaskPtr snptask,genetask;
  MKL_INT x_storage=VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT estimate;
  int errcode;
  hid_t file, dataset,filespace;
  hsize_t cordims[3]={snpsize,genesize,bsi};
  double* snpmeans,*snpvariances,*genemeans,*genevariances,*snp2r,*gene2r;
  

  C = (double *)mkl_malloc( snpsize*genesize*sizeof( double ), 64 );
  for (int i = 0; i < (snpsize*genesize); i++) {
    C[i] = 0.0;
  }
  readmatrix(Expmat,casetotal,genetotal,mycasestart,mycaseend,mygenestart,mygeneend,snpgenefilename,"genes");  
  readmatrix(Snpmat,casetotal,snptotal,mycasestart,mycaseend,mysnpstart,mysnpend,snpgenefilename,"snps");


  
  file = H5Fcreate(corfilename,H5F_ACC_RDWR,H5P_DEFAULT,H5P_DEFAULT);
  filespace = H5Screate_simple(3,cordims,NULL);
  dataset = H5Dcreate2(file,"cor",H5T_NATIVE_DOUBLE,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  

  //Bootstrapping loop
  //Initialize Bootstrap Matrices
  BootSNP= (double *) mkl_malloc(rowsize*snpsize*sizeof(double),64);
  BootExp= (double *) mkl_malloc(rowsize*genesize*sizeof(double),64);
  
  snpmeans = (double*) mkl_malloc(snpsize*sizeof(double),64);
  genemeans = (double*) mkl_malloc(genesize*sizeof(double),64);

  snpvariances = (double*) mkl_malloc(snpsize*sizeof(double),64);
  genevariances = (double*) mkl_malloc(genesize*sizeof(double),64);

  snp2r = (double*) mkl_malloc(snpsize*sizeof(double),64);
  gene2r = (double*) mkl_malloc(genesize*sizeof(double),64);

  
  
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,123);
  errcode = vsldSSNewTask(&snptask,&snpsize,&casesize,&x_storage,BootSNP,0,0);
  errcode = vsldSSEditTask(snptask,VSL_SS_ED_2R_MOM,snp2r);
  errcode = vsldSSEditTask(snptask,VSL_SS_ED_2C_MOM,snpvariances);
  errcode = vsldSSEditTask(snptask,VSL_SS_ED_MEAN,snpmeans);



  errcode = vsldSSNewTask(&genetask,&genesize,&casesize,&x_storage,BootExp,0,0);
  errcode = vsldSSEditTask(genetask,VSL_SS_ED_2R_MOM,gene2r);
  errcode = vsldSSEditTask(genetask,VSL_SS_ED_2C_MOM,genevariances);
  errcode = vsldSSEditTask(genetask,VSL_SS_ED_MEAN,genemeans);

  estimate = VSL_SS_MEAN|VSL_SS_2R_MOM|VSL_SS_2C_MOM;




  s_elapsed_boot=0;
  s_elapsed_mult=0;
  s_elapsed_write=0;
  for( int q=0; q<bsi; q++){
    cout<<"BSI:"<<q<<endl;
    s_initial_boot=dsecnd();
    DoBootstrap(BootSNP,BootExp,Snpmat,Expmat,casesize,snpsize,genesize,q,stream);
    vsldSSCompute(snptask,estimate,VSL_SS_METHOD_FAST);
    vdSqrt(snpsize,snpvariances,snpvariances,NULL);
    vdSqrt(genesize,genevariances,genevariances,NULL);
    vsldSSCompute(genetask,estimate,VSL_SS_METHOD_FAST);
    Scale(BootSNP,casesize,snpsize,snpmeans,snpvariances);
    Scale(BootExp,casesize,genesize,genemeans,genevariances);
    s_elapsed_boot+=(dsecnd()-s_initial_boot);
    s_initial_mult=dsecnd();
    //printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, snpsize, genesize, casesize, (1.0/(double) casesize), BootSNP, snpsize, BootExp, genesize, beta, C, genesize);
    s_elapsed_mult+=(dsecnd()-s_initial_mult);
    //printf ("\n Computations completed.\n\n");
    
    s_initial_write=dsecnd();
    WriteMat(C,mysnpstart,mysnpend,mygenestart,mygeneend,q,file,dataset);
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
  
  cordims[3]=2;

  file = H5Dcreate(annofilename,H5F_ACC_RDWR,H5P_DEFAULT,H5P_DEFAULT);
  filespace = H5Screate_simple(3,cordims,NULL);
  dataset = H5Dcreate2(file,"quantiles",H5T_NATIVE_DOUBLE,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  lowquantiles = (double*)mkl_malloc(snpsize*genesize*sizeof(double),64);
  highquantiles = (double*)mkl_malloc(snpsize*genesize*sizeof(double),64);


  for(int i=0; i<genesize; i++)
    {
      GetSlice(corfilename,mysnpstart,i,0,snpsize,1,bsi,bootmatrix);
      GetQuantile(bootmatrix,snpsize,bsi,tquantile);
      cblas_dcopy(snpsize,tquantile,2,&lowquantiles[index(0,i,genesize)],genesize);
      cblas_dcopy(snpsize,&tquantile[1],2,&highquantiles[index(1,i,genesize)],genesize);
    }
  WriteMat(lowquantiles,snpstart,snpend,genestart,geneend,0,dataset,filespace);
  WriteMat(highquantiles,snpstart,snpend,genestart,geneend,1,dataset,filespace);
  
  status = H5Dclose(dataset);
  status = H5Sclose(filespace);
  status = H5Fclose(file);  

  
      

  //Start Computing quantiles for correlation
  
  //printf (" Example completed. \n\n");
  return 0;
}
