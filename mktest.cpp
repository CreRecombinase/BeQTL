/* C source code is found in dgemm_example.c */

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
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
  int *snpchr,*snppos,*genechr,*genestart,*geneend;

  string *snpname,*genename;
  bool *cistrans;


  double mindist= 100000;
  hid_t filespace,dset_id;
  hsize_t* cordims;

  //file parameters
  char corfilename[]="/home/nwk2/mkl_test/testcor.nc";
  char snpfilename[]="/scratch/nwk2/hdf5/testsnp.h5";
  char genefilename[]="/scratch/nwk2/hdf5/testexp.h5";

  int casetotal = 177;
  int snptotal = 906598;
  int genetotal = 20501;

  int snpstart=0;
  int snpend=200;
  int casestart=0;
  int caseend=20;
  int genestart=0;
  int geneend=100;

  //Parameters for matrix multiplication using dgemm
  alpha = 1.0; beta = 0.0;  

  
  
  

  srand(time(0));
  

  
  C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
  for (i = 0; i < (m*n); i++) {
    C[i] = 0.0;
  }

  printf (" Intializing matrix data \n\n");
  
  readmatrix(Snpmat,casetotal,snptotal,casestart,caseend,snpstart,snpend,snpfilename,SNP);
  readmatrix(Expmat,casetotal,genetotal,casestsart,caseend,genestart,geneend,genefilename,GENE);

  //readanno(snpname,snpchr,snppos,NULL,m,"/home/nwk2/mkl_test/snpanno.txt");
  //readanno(genename,genechr,genestart,geneend,m,"/scratch/nwk2/mkl_test/expanno.txt");
  /* Let's establish cis and trans relationships now*/
  /*
  //  cistrans = (bool *)mkl_malloc(n*m*sizeof(bool),64);
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
