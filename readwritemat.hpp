#ifndef READWRITEMAT_HPP
#define READWRITEMAT_HPP
#include "hdf5.h"
#include "hdf5_hl.h"
#include<iostream>
#include <fstream>
#include <cstring>
#include "mktest.hpp"
#include "mkl.h"
#include <sys/stat.h>
#include <map>


struct cmp_str
{
  bool operator()(char const *a, char const *b)
  {
    return(strcmp(a,b)<0);
  }
};

//Return Distance instead of a boolean
int CisDist(const int* snparray,const int* genearray,int snp,int gene)
{
  int snpchr = snparray[index(snp,0,2)];
  int snppos = snparray[index(snp,1,2)];
  int genechr = genearray[index(gene,0,3)];
  int genestart = genearray[index(gene,1,3)];
  int geneend = genearray[index(gene,2,3)];
  if(snpchr!=genechr)
    {
      return(-1);
    }
  return(min(abs(snppos-genestart),abs(snppos-geneend)));
}
/*  
  return( (snpchr==genechr)&&((abs(snppos-genestart)<cisdist)||(abs(snppos-geneend)<cisdist)));
  
}
*/


void PrintMat(const double *matrix,int r,int c, const char* filename)
{
  ofstream file;
  file.open(filename,ofstream::out|ofstream::app);
  for(int i=0;i<c;i++)
    {
      for(int j=0;j<r;j++)
	{
	  file<<matrix[index(i,j,c)];
	  if(j!=r)
	    file<<"\t";
	}
      file<<"\n";
    }
  file.close();
  //    cout<<"____________________"<<endl;
}
void PrintMat(const int *matrix,int r,int c, const char* filename)
{
  ofstream file;
  file.open(filename,ofstream::out);
  for(int i=0;i<r;i++)
    {
      for(int j=0;j<c;j++)
	{
	  file<<matrix[index(i,j,c)];
	  if(j!=c)
	    file<<"\t";
	}
      file<<"\n";
    }
  file.close();
  //    cout<<"____________________"<<endl;
}

void PrintMat(const double *matrix,int r,int c)
//Function for printing 2D, row order array
//*matrix   pointer to head of matrix
//r         number of rows
//c         number of columns
{
  for(int i=0;i<r;i++)
    {
      for(int j=0;j<c;j++)
	{
	  cout<<matrix[index(i,j,c)];
	  if(j!=c)
	    cout<<"\t";
	}
      cout<<"\n";
    }
    cout<<"____________________"<<endl;
}

void PrintMat(const int *matrix,int r,int c)
//Function for printing 2D, row order array
//*matrix   pointer to head of matrix
//r         number of rows
//c         number of columns
{
  for(int i=0;i<r;i++)
    {
      for(int j=0;j<c;j++)
	{
	  cout<<matrix[index(i,j,c)];
	  if(j!=c)
	    cout<<"\t";
	}
      cout<<"\n";
    }
    cout<<"____________________"<<endl;
}



herr_t ReadMatrix(int *&matrix,int rstart,int rsize,int cstart,int csize,const char* filename, const char* snpgene)
{

  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;
  hsize_t dimsout[2]={rsize,csize};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={rsize,csize};
    
  file  = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset = H5Dopen2(file,snpgene,H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);

  memspace = H5Screate_simple(2,dimsm,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);

  status = H5Dread(dataset,H5T_NATIVE_INT,memspace,dataspace,H5P_DEFAULT,matrix);

  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);

  return(status);
}

herr_t ReadMatrix(char** &matrix,int rstart,int rsize,const char* filename, const char* snpgene)
{

 
  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;
 
  hsize_t dimsout[1]={rsize};
  hsize_t offset[1]={rstart};
  hsize_t offset_out[1] ={0};
  hsize_t dimsm[1] ={rsize};
  hid_t fstrtype;
   


  

 

  
  file  = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  dataset = H5Dopen2(file,snpgene,H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);

  memspace = H5Screate_simple(1,dimsm,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);

  status = H5Dread(dataset,fstrtype,memspace,dataspace,H5P_DEFAULT,matrix);

  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);

  return(status);
}

herr_t ReadMatrix(double *&matrix,const int rstart,const int rsize, const int cstart, const int csize,const char* filename,const char* snpgene)
{

  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;

  hsize_t dimsout[2]={rsize,csize};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={rsize,csize};

   


  


  
  file  = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset = H5Dopen2(file,snpgene,H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);

  memspace = H5Screate_simple(2,dimsm,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);

  status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,matrix);

  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);

  return(status);




}


herr_t readparam(hid_t file_id, const char* attributename,char*&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  herr_t err;

  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  fstrtype= H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);
  att=H5Aopen_name(root,attributename);
  err=H5Aread(att,fstrtype,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
  return(err);
  
}
	  
herr_t readparam(hid_t file_id,const char* attributename,int&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  herr_t err;
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  err=H5Aread(att,H5T_NATIVE_INT,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
  return(err);
}

herr_t readparam(hid_t file_id,const char* attributename,double&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  herr_t err;
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  err=H5Aread(att,H5T_NATIVE_DOUBLE,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
  return(err);
}

map<char*,int,cmp_str> ReadAnno(const char* annofile,const int totalnum, const char* snpgene)
{
  char **names;

  map<char*, int, cmp_str> mmap;
  names=(char**) mkl_malloc(totalnum*sizeof(char*),64);
  ReadMatrix(names,0,totalnum,annofile,snpgene);

  for(int i=0; i<totalnum; i++)
    {
      mmap[names[i]]=i;
    }
  mkl_free(names);
  return(mmap);
}



void CisTransOut (const double *quantilearray,const int snpstart,const int snpsize,const int genestart,const int genesize, const int casesize,const int zsize,const char *annofile,double t_thresh,const char* eqtlfile,int mychunk,const char* progfile,map<char*,int,cmp_str>snpanno, map<char*,int,cmp_str> geneanno, const int *snparray, const int *genearray, const char* snpgenefile)
  
{
  
  struct stat buffer;
  int status;
  int cisdist;
  FILE* fp;
  ofstream eqtl,prog;
  double cor1,cor2,cor3;
  double cordenom;

  char**dsnpnames,**dgenenames;

  
  cordenom = sqrt(casesize-2);
  dsnpnames= (char**)mkl_malloc(snpsize*sizeof(char*),64);
  dgenenames=(char**)mkl_malloc(genesize*sizeof(char*),64);
  ReadMatrix(dsnpnames,snpstart,snpsize,snpgenefile,"snpnames");

  ReadMatrix(dgenenames,genestart,genesize,snpgenefile,"genenames");
  
  
  prog.open(progfile,ofstream::out|ofstream::app);
  prog<<mychunk<<"\t"<<snpsize*genesize<<endl;
  eqtl.open(eqtlfile,ofstream::out|ofstream::app);
  if(mychunk==0){
    eqtl<<"SNP\tGene\tt-statLow\tt-statMed\tt-statHigh\tDistance"<<endl;
  }
  
  for(int i=0; i<snpsize; i++)
    {
      for(int j=0; j<genesize; j++)
	{
	  cor1=quantilearray[multindex(i,j,0,snpsize,genesize,zsize)];
	  cor1=cordenom*(cor1/sqrt(1-(cor1*cor1)));
	  cor2=quantilearray[multindex(i,j,1,snpsize,genesize,zsize)];
	  cor2=cordenom*(cor2/sqrt(1-(cor2*cor2)));
	  cor3=quantilearray[multindex(i,j,2,snpsize,genesize,zsize)];
	  cor3=cordenom*(cor3/sqrt(1-(cor3*cor3)));
	  //	  cout<<snpnames[i]<<"-"<<genenames[j]<<": "<<fabs(worstcor)<<" "<<t_thresh<<endl;
	  
	  if(fabs(cor2)>t_thresh)
	    {
	      cisdist = CisDist(snparray,genearray,snpanno[dsnpnames[i]],geneanno[dgenenames[j]]);
	      eqtl<<dsnpnames[i]<<"\t"<<dgenenames[j]<<"\t"<<cor1<<"\t"<<cor2<<"\t"<<cor3<<"\t"<<cisdist<<endl;
	      
	    }
	}
    }
	    
  mkl_free(dsnpnames);
  mkl_free(dgenenames);
  prog.close();
  eqtl.close();
}


/*
void WriteQuantiles(const char* filename,const double* quantarray,const size_t snpstart,const size_t snpsize,const size_t genestart,const size_t genesize,const int mychunk,const char* lockfile,const char* snpgenefilename)
{
  struct stat buffer;
  FILE* fp;
  int snptotal,genetotal,casesize;
  double cor,cordenom;
  double* B;

  hid_t file_id;
  hid_t cordataspace,cordataset,cormemspace;
  herr_t status;
  hsize_t cor_max_size[3];
  hsize_t cor_sl_size[3]={snpsize,genesize,1};
  hsize_t cor_sl_offset[3]={snpstart,genestart,0};

  file_id = H5Fopen(snpgenefilename,H5F_ACC_RDONLY,H5P_DEFAULT);
  readparam(file_id,"snptotal",snptotal);
  readparam(file_id,"genetotal",genetotal);
  readparam(file_id,"casetotal",casesize);
  H5Fclose(file_id);
  
  cor_max_size[0]=snptotal;
  cor_max_size[1]=genetotal;
  cor_max_size[2]=2;




  
      

  cordenom = sqrt(casesize-2);
  B=(double*)mkl_malloc(2*snpsize*genesize*sizeof(double),64);

  vdpowx(2*snpsize*genesize,quantarray,2,B);
  dscal(2*snpsize*genesize,quantarray,cordenom,1);
  dscal(2*snpsize*genesize,(-1.0),B,1);
  for(int i=0; i<2*snpsize*genesize;i++)
    {
      ++(B[i]);
    }
  vdsqrt(2*snpsize*genesize,B,B);
  vddiv(2*snpsize*genesize,quantarray,B,quantarray);
  
  /*
  for(size_t i=0; i<snpsize; i++)
    {
      for(size_t j=0; j<genesize; j++)
	{
	  cor = quantarray[multindex(i,j,0,snpsize,genesize,2)];
	  corarray[index(i,j,genesize)]=cordenom*(cor/sqrt(1-(cor*cor)));
	}
    }
  while(stat(lockfile,&buffer)==0)
    {
      usleep(999999);
      cout<<lockfile<<" exists!"<<endl;
    }
  fp = fopen(lockfile,"ab+");
  if(mychunk==0)
    {
      file_id=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      cordataspace=H5Screate_simple(3,cor_max_size,NULL);
      cordataset=H5Dcreate(file_id,"corquantiles",H5T_NATIVE_DOUBLE,cordataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      H5Sclose(cordataspace);
      H5Dclose(cordataset);
      H5Fclose(file_id);
    }
  
  file_id = H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
  cormemspace = H5Screate_simple(3,cor_sl_size,NULL);
  cordataset = H5Dopen2(file_id,"corquantiles",H5P_DEFAULT);
  cordataspace=H5Dget_space(cordataset);
  H5Sselect_hyperslab(cordataspace,H5S_SELECT_SET,cor_sl_offset,NULL,cor_sl_size,NULL);
  H5Dwrite(cordataset,H5T_NATIVE_DOUBLE,cormemspace,cordataspace,H5P_DEFAULT,corarray);
  
  H5Dclose(cordataset);
  H5Sclose(cormemspace);
  H5Sclose(cordataspace);
  H5Fclose(file_id);
  fclose(fp);
  remove(lockfile);
  
  
  for(size_t i=0; i<snpsize; i++)
    {
      for(size_t j=0; j<genesize; j++)
	{
	  cor = quantarray[multindex(i,j,1,snpsize,genesize,2)];
	  corarray[index(i,j,genesize)]=cordenom*(cor/sqrt(1-(cor*cor)));
	}
    }

  while(stat(lockfile,&buffer)==0)
    {
      usleep(999999);
      cout<<lockfile<<" exists!"<<endl;
    }
  fp = fopen(lockfile,"ab+");

  cor_sl_offset[2]=1;
  file_id = H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
  cormemspace = H5Screate_simple(3,cor_sl_size,NULL);
  cordataset = H5Dopen2(file_id,"corquantiles",H5P_DEFAULT);
  cordataspace=H5Dget_space(cordataset);
  H5Sselect_hyperslab(cordataspace,H5S_SELECT_SET,cor_sl_offset,NULL,cor_sl_size,NULL);
  H5Dwrite(cordataset,H5T_NATIVE_DOUBLE,cormemspace,cordataspace,H5P_DEFAULT,corarray);
  
  H5Dclose(cordataset);
  H5Sclose(cordataspace);
  H5Sclose(cormemspace);
  H5Fclose(file_id);
  fclose(fp);
  remove(lockfile);

  mkl_free(corarray);
}
*/
//void 


#endif
