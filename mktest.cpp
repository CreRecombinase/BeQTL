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
  if(argc<6)
    {
      cerr<<"usage: "<<argv[0]<<" startchunk endchunk snpgenefilename readfilelock writefilelock progfile"<<endl;
      return(-1);
    }
  
  //Large matrices for main loop
  double *Snpmat=NULL, *Genemat=NULL, *C=NULL;
  double *BootSNP=NULL, *BootGene=NULL;

 
  //Matrices and temporary vectors for reporting and storing quantiles
  double *tquantile;
  

  //file parameters

  char *snpgenefilename,*annofilename,*cisfilename,*transfilename,*readfilelock,*writefilelock,*progfile;


  // Platform data/case numbers
  int casetotal;
  int snptotal;
  int genetotal;

  //variables for profiling
  double s_initial_boot,s_elapsed_boot;
  double s_initial_mult, s_elapsed_mult;
  double s_initial_quant,s_elapsed_quant;
  double s_initial_write,s_elapsed_write;
  double s_initial_memset,s_elapsed_memset;
  double s_initial_read,s_elapsed_read;

  int mychunk,lastchunk,mysnpchunk,mygenechunk,chunkstart;
  int snpchunks;
  int genechunks;

  int snpsize;
  int genesize;

  int casesize;
  int snpstart;
  int casestart=0;
  int genestart;
  int cisdist;
  int bsi;
  int alignsize=64;
  
  double alpha,beta,t_thresh,mean=0;
  int  *bootrows;
  size_t C_dim=0,new_c_dim=0;
  
  chunkstart = atoi(argv[1]);
  lastchunk = atoi(argv[2]);
  snpgenefilename=argv[3];
  readfilelock=argv[4];
  writefilelock=argv[5];
  progfile=argv[6];
  
  hid_t file_id;
  file_id = H5Fopen(snpgenefilename,H5F_ACC_RDONLY,H5P_DEFAULT);

  readparam(file_id,"annofile",annofilename);
  readparam(file_id,"cisfilename",cisfilename);
  readparam(file_id,"transfilename",transfilename);
  readparam(file_id,"snpchunks",snpchunks);
  readparam(file_id,"genechunks",genechunks);
  readparam(file_id,"snptotal",snptotal);
  readparam(file_id,"genetotal",genetotal);
  readparam(file_id,"casetotal",casesize);
  readparam(file_id,"snpsize",snpsize);
  readparam(file_id,"genesize",genesize);
  casetotal=casesize;
  readparam(file_id,"bsi",bsi);

  readparam(file_id,"t_thresh",t_thresh);
  readparam(file_id,"cisdist",cisdist);

  H5Fclose(file_id);


  cout<<"parameterization successful!"<<endl;
		


  alpha = 1.0/((double)casesize);
  beta = 0.0;  
  s_elapsed_read=0;
  s_elapsed_boot=0;
  s_elapsed_mult=0;
  s_elapsed_quant=0;
  s_elapsed_write=0;
  
  for(mychunk =chunkstart; mychunk<=lastchunk; mychunk++)
    {
      if(mychunk>snpchunks*genechunks)
	{
	  cout<<"eQTL analysis complete!"<<endl;
	  return(0);
	}

      //parameters set from mychunk
      mysnpchunk  = mychunk/(genechunks);
      mygenechunk = mychunk%(genechunks);

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
      new_c_dim = (size_t)snpsize*(size_t)genesize*(size_t)bsi*sizeof(double);


  
      cout<<"Initializing large corarray of size:"<<snpsize<<"*"<<genesize<<"*"<<bsi<<"*"<<sizeof(double)<<"="<<new_c_dim<<endl;
      if(new_c_dim>C_dim)
	{
	  if(C!=NULL)
	    {
	      mkl_free(C);
	      mkl_free(BootSNP);
	      mkl_free(BootGene);
	    }
	  C_dim = (size_t)snpsize*(size_t)genesize*(size_t)bsi*sizeof(double);
	  C = (double *)mkl_malloc( C_dim,alignsize);
	  BootSNP= (double *) mkl_malloc(casesize*snpsize*sizeof(double),alignsize);
	  BootGene= (double *) mkl_malloc(casesize*genesize*sizeof(double),alignsize);
	}

      if(C==NULL)
	{

	  cout<<"Bad things happened!"<<endl;
	}
      C[(C_dim/sizeof(double))-1]=0;

    
      cout<<"SNPrange= "<<snpstart<<"-"<<snpstart+snpsize<<endl;
      cout<<"GENErange= "<<genestart<<"-"<<genestart+genesize<<endl;
      
      s_initial_read=dsecnd();
      cout<<"Reading in gene matrix!"<<endl;
      ReadMatrix(Genemat,casetotal,genetotal,casestart,casesize,genestart,genesize,snpgenefilename,"genes",readfilelock);
      cblas_daxpby(casesize*genesize,100,Genemat,1,0,Genemat,1);
	  
      cout<<"Reading in SNP matrix!"<<endl;
      ReadMatrix(Snpmat,casetotal,snptotal,casestart,casesize,snpstart,snpsize,snpgenefilename,"snps",readfilelock);
   
  
      s_elapsed_read += dsecnd()-s_initial_read;
      
      //Bootstrapping loop
      //Initialize Bootstrap Matrices
     
      
      cout<<"Making bootstrap"<<endl;
      if(mychunk==chunkstart)
	{
	  MakeBootRows(bootrows,casesize,bsi);
	}

      alpha = 1/((double) casesize-1);

      cout<<"Starting Bootstrap"<<endl;
      //  s_elapsed_write=0;
      for( size_t q=0; q<bsi; q++){
	int percent = (((double) q)/((double) bsi))*100;
	printProgBar(percent);
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
	
	//	cout<<"scaling done; about to mult"<<endl;
	s_initial_mult=dsecnd();
	
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, snpsize, genesize, casesize,1, BootSNP, snpsize, BootGene, genesize, beta, &C[ci], genesize);
	s_elapsed_mult+=(dsecnd()-s_initial_mult);
	

	//PrintMat(&C[(size_t)(snpsize*genesize*q)],snpsize,genesize);
	//	PrintMat(C,snpsize,genesize);
      }
      
      //            cout<<"Deallocating memory"<<endl;
      
      
      mkl_free(Snpmat);
      mkl_free(Genemat);
      if(mychunk==lastchunk)
	{
	  mkl_free(BootSNP);
	  mkl_free(BootGene);
	}
 
      //Initializing quantile matrix
      //    cout<<"Allocating quantile matrix"<<endl;
      tquantile = (double*)mkl_malloc(2*(size_t)snpsize*(size_t)genesize*sizeof(double),alignsize);

      //    cout<<"Computing quantiles"<<endl;
      s_initial_quant=dsecnd();
      GetQuantile(C,snpsize,genesize,bsi,tquantile);
      s_elapsed_quant+=dsecnd()-s_initial_quant;

      if(mychunk==lastchunk)
	{
	  mkl_free(C);
	}

      cout<<"Writing quantiles to file: "<<cisfilename<<endl;
      s_initial_write=dsecnd();
      CisTransOut(tquantile,snpstart,snpsize,genestart,genesize,casesize,2,annofilename,t_thresh,cisfilename,transfilename,cisdist,writefilelock,mychunk,progfile);
      s_elapsed_write+=dsecnd()-s_initial_write;
      mkl_free(tquantile);
      cout<<"Read time: "<<s_elapsed_read<<endl;
      cout<<"Bootstrap time: "<<s_elapsed_boot<<endl;
      cout<<"Mult time: "<<s_elapsed_mult<<endl;
      cout<<"Quantile time: "<<s_elapsed_quant<<endl;
      cout<<"Write time: "<<s_elapsed_write<<endl;
    }

  return (0);
}
