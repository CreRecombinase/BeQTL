#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl.h"
#include "mktest.hpp"
#include "readwritemat.hpp"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <sstream>
#include "strideiter.hpp"



using namespace std;

int main(int argc, char* argv[])
{
  if(argc<5)
    {
      cerr<<"usage: "<<argv[0]<<" startchunk numchunks snpgenefilename threads progbar"<<endl;
      return(-1);
    }
  
  //Large matrices for main loop
  double *Snpmat=NULL, *Genemat=NULL, *C=NULL;
  double *BootSNP=NULL, *BootGene=NULL;

 
  //Matrices and temporary vectors for reporting and storing quantiles
  double *tquantile;
  int* snpannoarray,*geneannoarray;

  

  //file parameters
  char *snpgenefilename,*annofilename,*teqtlfilename,*tprogfilename;
  


  // Platform data/case numbers
  int progbar;
  int threads;
  int snptotal,genetotal;
  int snpabstotal,geneabstotal;

  //variables for profiling
  double s_initial_boot,s_elapsed_boot;
  double s_initial_mult, s_elapsed_mult;
  double s_initial_quant,s_elapsed_quant;
  double s_initial_quant2,s_elapsed_quant2;
  double s_initial_write,s_elapsed_write;
  double s_initial_memset,s_elapsed_memset;
  double s_initial_read,s_elapsed_read;

  int mychunk,lastchunk,mysnpchunk,mygenechunk,chunkstart;
  int numchunks;
  int snpchunks;
  int genechunks;

  int osnpsize;
  int ogenesize;
  size_t snpsize;
  size_t genesize;

  int casesize;
  int snpstart;
  int casestart=0;
  int genestart;
  int bsi;
  int alignsize=64;
  int quantilechunks;
  
  double alpha,beta,t_thresh;
  int  *bootrows;
  size_t C_dim=0,new_c_dim=0;
  
  chunkstart = atoi(argv[1]);
  numchunks = atoi(argv[2]);
  lastchunk=chunkstart+numchunks;
  snpgenefilename=argv[3];
  threads=atoi(argv[4]);
  threads=threads*2;
  progbar=atoi(argv[5]);
  
  hid_t file_id;
  herr_t err=0;
  file_id = H5Fopen(snpgenefilename,H5F_ACC_RDONLY,H5P_DEFAULT);

  err=err|readparam(file_id,"annofile",annofilename);
  err=err|readparam(file_id,"eqtlfilename",teqtlfilename);
  err=err|readparam(file_id,"progfile",tprogfilename);
  err=err|readparam(file_id,"snpchunks",snpchunks);
  err=err|readparam(file_id,"genechunks",genechunks);
  err=err|readparam(file_id,"snptotal",snptotal);
  err=err|readparam(file_id,"genetotal",genetotal);
  err=err|readparam(file_id,"casetotal",casesize);
  err=err|readparam(file_id,"snpsize",osnpsize);
  snpsize=osnpsize;
  err=err|readparam(file_id,"genesize",ogenesize);
  genesize=ogenesize;
  err=err|readparam(file_id,"bsi",bsi);
  err=err|readparam(file_id,"snpabstotal",snpabstotal);
  err=err|readparam(file_id,"geneabstotal",geneabstotal);

  err=err|readparam(file_id,"t_thresh",t_thresh);
  if(err<0){
    cerr<<"Parameter not successfully read!: "<<err<<endl;
    return(err);
  }

  H5Fclose(file_id);



  cout<<"parameterization successful!"<<endl;
		


  alpha = 1.0/((double)casesize);

  beta = 0.0;  
  s_elapsed_read=0;
  s_elapsed_boot=0;
  s_elapsed_mult=0;
  s_elapsed_quant=0;
  s_elapsed_write=0;

  //Set number of threads to be number of slots on LSF
  mkl_set_num_threads(threads);

  cout<<"Initializing correlation array of size: "<<snpsize<<"*"<<genesize<<"*"<<bsi<<"*"<<sizeof(double)<<"="<<new_c_dim<<endl; 
  C_dim = (size_t)snpsize*(size_t)genesize*(size_t)bsi*sizeof(double);
  C = (double *)mkl_malloc( C_dim,alignsize);
  Snpmat = (double*)mkl_malloc(casesize*snpsize*sizeof(double),alignsize);
  Genemat = (double*)mkl_malloc(casesize*snpsize*sizeof(double),alignsize);
  BootSNP= (double *) mkl_malloc(casesize*snpsize*sizeof(double),alignsize);
  BootGene= (double *) mkl_malloc(casesize*genesize*sizeof(double),alignsize);  
  tquantile = (double*)mkl_malloc(3*(size_t)snpsize*(size_t)genesize*sizeof(double),alignsize);
  bootrows = (int*) mkl_malloc((size_t)casesize*(size_t)bsi*sizeof(int),alignsize);
  snpannoarray= (int*) mkl_malloc(snpabstotal*2*sizeof(int),alignsize);
  geneannoarray= (int*) mkl_malloc(geneabstotal*3*sizeof(int),alignsize);


  ReadMatrix(snpannoarray,0,snpabstotal,0,2,annofilename,"snparray");
  ReadMatrix(geneannoarray,0,geneabstotal,0,3,annofilename,"genearray");

  
  map<char*,int,cmp_str> snpanno=ReadAnno(annofilename,2,"snpname");
  map<char*,int,cmp_str> geneanno=ReadAnno(annofilename,3,"genename");


  if(C==NULL)
    {
      
      cout<<"Not enough memory for corarray!"<<endl;
      return(-1);
    }
  C[(C_dim/sizeof(double))-1]=0;
  stringstream sss;
  sss<<chunkstart;
  string chunkstr=sss.str();
  string progfilename=std::string(tprogfilename)+std::string(chunkstr)+".txt";
  string eqtlfilename=std::string(teqtlfilename)+std::string(chunkstr)+".txt";
  for(mychunk =chunkstart; mychunk<lastchunk; mychunk++)
    {


     
      if(mychunk>snpchunks*genechunks)
	{
	  cout<<"eQTL analysis complete!"<<endl;
	  return(0);
	}

      //parameters set from mychunk
      mysnpchunk  = mychunk/(genechunks);
      mygenechunk = mychunk%(genechunks);
      snpsize=osnpsize;
      genesize=ogenesize;
      snpstart =mysnpchunk*snpsize;
      if(snpsize+snpstart>snptotal)
	{
	  snpsize=snptotal-snpstart;
	}
      genestart = mygenechunk*genesize;
      if(genesize+genestart>genetotal)
	{
	  genesize = genetotal-genestart;
	}
      cout<<"iteration number: "<<mychunk<<" SNPchunk: "<<mysnpchunk<<" Genechunk: "<<mygenechunk<<endl;
      //Size of result array

    
      cout<<"SNPrange= "<<snpstart<<"-"<<snpstart+snpsize<<endl;
      cout<<"GENErange= "<<genestart<<"-"<<genestart+genesize<<endl;
      
      s_initial_read=dsecnd();
      cout<<"Reading in gene matrix!"<<endl;
      ReadMatrix(Genemat,casestart,casesize,genestart,genesize,snpgenefilename,"genes");
      cblas_dscal(casesize*genesize,100,Genemat,1);
	  
      cout<<"Reading in SNP matrix!"<<endl;
      ReadMatrix(Snpmat,casestart,casesize,snpstart,snpsize,snpgenefilename,"snps");
      s_elapsed_read += dsecnd()-s_initial_read;
      
      //Bootstrapping loop
      //Initialize Bootstrap Matrices
     
      
      if(mychunk==chunkstart)
	{
	  MakeBootRows(bootrows,casesize,bsi);
	}



      cout<<"Starting Bootstrap"<<endl;
      for( size_t q=0; q<bsi; q++){
	int percent = (((double) q)/((double) bsi))*100;
	if(progbar>0)
	  {
	    printProgBar(percent);
	  }
	//	cout<<q<<endl;
	s_initial_boot=dsecnd();
	DoBootstrap(BootSNP,Snpmat,casesize,snpsize,&bootrows[index(q,0,casesize)]);
	DoBootstrap(BootGene,Genemat,casesize,genesize,&bootrows[index(q,0,casesize)]);
	DoScale(BootSNP,casesize,snpsize);
	DoScale(BootGene,casesize,genesize);
	cblas_daxpby(casesize*genesize,alpha,BootGene,1,0,BootGene,1);
	s_elapsed_boot+=(dsecnd()-s_initial_boot);
	size_t ci = snpsize*genesize*q;
	//	cout<<ci<<"\t"<<C_dim<<endl;
	

	s_initial_mult=dsecnd();
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, snpsize, genesize, casesize,1, BootSNP, snpsize, BootGene, genesize, beta, &C[ci], genesize);
	s_elapsed_mult+=(dsecnd()-s_initial_mult);
	
      }
      

 
      //Initializing quantile matrix
      //    cout<<"Allocating quantile matrix"<<endl;


      cout<<"Computing quantiles"<<endl;


      s_initial_quant=dsecnd();
      GetQuantile(C,snpsize,genesize,bsi,tquantile);
      //cout<<"Writing quantiles to "<<corfilename<<endl;
      //WriteQuantiles(corfilename,tquantile,snpstart,snpsize,genestart,genesize,mychunk,writefilelock,snpgenefilename);
      s_elapsed_quant+=dsecnd()-s_initial_quant;
      cout<<"Writing quantiles to file: "<<eqtlfilename<<endl;
      s_initial_write=dsecnd();
      CisTransOut(tquantile,snpstart,snpsize,genestart,genesize,casesize,3,annofilename,t_thresh,eqtlfilename.c_str(),mychunk,progfilename.c_str(),snpanno,geneanno,snpannoarray,geneannoarray,snpgenefilename);
      s_elapsed_write+=dsecnd()-s_initial_write;
      cout<<"Read time: "<<s_elapsed_read<<endl;
      cout<<"Bootstrap time: "<<s_elapsed_boot<<endl;
      cout<<"Mult time: "<<s_elapsed_mult<<endl;
      cout<<"Quantile time: "<<s_elapsed_quant<<endl;
      cout<<"Write time: "<<s_elapsed_write<<endl;


      s_elapsed_read=0;
      s_elapsed_boot=0;
      s_elapsed_mult=0;
      s_elapsed_quant=0;
      s_elapsed_write=0;
    }
      mkl_free(Snpmat);
      mkl_free(bootrows);
      mkl_free(tquantile);
      mkl_free(Genemat);
      mkl_free(C);
      mkl_free(BootSNP);
      mkl_free(BootGene);


  return (0);
}
