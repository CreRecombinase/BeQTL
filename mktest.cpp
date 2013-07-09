/* C source code is found in dgemm_example.c */

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "mkl.h"
#include "mktest.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <time.h>





using namespace std;

int main()
{
  double *Snpmat=NULL, *Expmat=NULL, *C=NULL;
  double *BootSNP=NULL, *BootExp=NULL;
  int *snpchr,*snppos,*genechr,*genestart,*geneend;

  string *snpname,*genename;
  bool *cistrans;
  
  int m, n, p, i, j;
  int bsi;
  double alpha, beta;
  double colmean=0;
  double colsd=0;
  double mindist= 100000;
  ifstream snpfile;
  ifstream expfile;
  string line,ts;

  hid_t file_id,dset_id;
  hid_t filespace,memspace;
  hsize_t dimsf[3];

  herr_t status;

  char filename[]="/home/nwk2/mkl_test/testcor.nc";

  

  srand(time(0));
 
  
  m = 20, p = 81, n = 10,bsi=4;
  dimsf[0]=m;
  dimsf[1]=n;
  dimsf[2]=bsi;
  file_id = H5Fcreate (filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  filespace = H5Screate_simple(3,dimsf,NULL);
  
  dset_id = H5Dcreate2(file_id,"DoubleArray",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Sclose(filespace);


  printf (" Initializing data for matrix multiplication C=Snpmat*Expmat for matrix \n"
	  " Snpmat(%ix%i) and matrix Expmat(%ix%i)\n\n", m, p, p, n);
  alpha = 1.0; beta = 0.0;
  
  printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
	  " performance \n\n");
  
  //  cistrans = (bool *)mkl_malloc(n*m*sizeof(bool),64);

  
  C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
  for (i = 0; i < (m*n); i++) {
    C[i] = 0.0;
  }

  printf (" Intializing matrix data \n\n");
  
  readmatrix(Snpmat,p,m,"/scratch/nwk2/hdf5/testsnp.txt");
  readmatrix(Expmat,p,n,"/scratch/nwk2/hdf5/testexp.txt");

  //readanno(snpname,snpchr,snppos,NULL,m,"/home/nwk2/mkl_test/snpanno.txt");
  //readanno(genename,genechr,genestart,geneend,m,"/scratch/nwk2/mkl_test/expanno.txt");
  /* Let's establish cis and trans relationships now*/
  /*
    for(i=0; i<m; i++){
    for(j=0; j<n; j++){
    cistrans[index(i,j,m)]=(snpchr[i]==genechr[j])&&((abs(snppos[i]-genestart[j])<mindist)||(abs(snppos[i]-geneend[j])<mindist));
    if(cistrans[index(i,j,m)])
    {
    cout<<"Cis found!"<<snpname[i]<<" "<<genename[j]<<endl;
    }
    }
    }
  */
  

  for( int q=0; q<bsi; q++){
    bootstrap(BootSNP,BootExp,Snpmat,Expmat,p,m,n,q);
    scale(BootSNP,p,m);
    scale(BootExp,p,n);
    
    printf (" Computing matrix product using Intel® MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 
		m, n, p, (1.0/(double) p), BootSNP, m, BootExp, n, beta, C, n);
    printf ("\n Computations completed.\n\n");
    
    writemat(C,m,n,filespace,dset_id,q);
  }
  printf ("\n Top left corner of matrix C: \n");
  for (i=0; i<10; i++) {
    for (j=0; j<10; j++) {
      printf ("%0.6f\t", C[index(i,j,n)]/p);
      
    }
    printf ("\n");
  }
  
  
  printf ("\n Deallocating memory \n\n");
  mkl_free(Snpmat);
  mkl_free(Expmat);
  mkl_free(C);
  
  printf (" Example completed. \n\n");
  return 0;
}
