#ifndef MKTEST_HPP
#define MKTEST_HPP
#include <string>
#include <fstream>
#include "hdf5.h"
#include <math.h>
#include "mkl.h"
#include "hdf5_hl.h"
#include <sys/file.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>



bool debug=false;
using namespace std;
enum snp_exp {SNP,GENE};

size_t multindex(size_t row,size_t column,size_t zindex,size_t nrows,size_t ncols,size_t zsize)
{
  return(row*(ncols*zsize)+(zsize*column)+zindex);
}

size_t index (size_t row, size_t column, size_t nocols)
//Function for indexing 2 dimensional, row major order array
//row     row index
//column  column index
//nocols  number of columns
{
 
  return((nocols*row)+column);
}

void PrintMat(const double *matrix,int r,int c, const char* filename)
{
  ofstream file;
  file.open(filename,ofstream::out);
  for(int i=0;i<r;i++)
    {
      for(int j=0;j<c;j++)
	{
	  file<<matrix[index(i,j,c)]<<"\t";
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
	  cout<<matrix[index(i,j,c)]<<"\t";
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
	  cout<<matrix[index(i,j,c)]<<"\t";
	}
      cout<<"\n";
    }
    cout<<"____________________"<<endl;
}


void DoScale(double *&matrix,int rows,int cols)
/*Function to perform scaling on the columns of the matrix such that the rows have a mean of 0 and unit variance
Scaling is performed by subtracting the column means and then dividing by the square root of the variance (the standard deviation)

matrix     pointer to head of array
rows       number of rows
cols       number of columns
 */
{
  VSLSSTaskPtr task;
  MKL_INT x_storage=VSL_SS_MATRIX_STORAGE_COLS;
  unsigned MKL_INT estimate;
  double *means,*variances,*secondraw;
  int errcode;

  means = (double*) mkl_malloc(cols*sizeof(double),64);
  variances = (double*)mkl_malloc(cols*sizeof(double),64);
  secondraw = (double*)mkl_malloc(cols*sizeof(double),64);
  

  errcode = vsldSSNewTask(&task,&cols,&rows,&x_storage,matrix,0,0);
  errcode = vsldSSEditTask(task,VSL_SS_ED_2R_MOM,secondraw);
  errcode = vsldSSEditTask(task,VSL_SS_ED_2C_MOM,variances);
  errcode = vsldSSEditTask(task,VSL_SS_ED_MEAN,means);

  estimate = VSL_SS_MEAN|VSL_SS_2R_MOM|VSL_SS_2C_MOM;
  errcode = vsldSSCompute(task,estimate,VSL_SS_METHOD_FAST);
  errcode = vslSSDeleteTask(&task);
  
  vdSqrt(cols,variances,variances);
  for(int i=0; i<rows; i++)
    {
  vdSub(cols,&matrix[index(i,0,cols)],means,&matrix[index(i,0,cols)]);
  vdDiv(cols,&matrix[index(i,0,cols)],variances,&matrix[index(i,0,cols)]);
    }
  mkl_free(means);
  mkl_free(variances);
  mkl_free(secondraw);

}

void MakeBootRows(int *&bootrows, int rowsize,int bstotal)
//Function to geneerate 2d array of bootstrap rows
{
  VSLStreamStatePtr stream;
  int errcode;
  
  bootrows = (int*) mkl_malloc((size_t)rowsize*(size_t)bstotal*sizeof(int),64);
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,123);
  errcode = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,rowsize*bstotal,bootrows,0,rowsize);
  errcode = vslDeleteStream(&stream);

  

}

void DoBootstrap(double *&bootmatrix,const double* omatrix,int rowsize, int colsize,int* bootrows)
/*Function that takes original array and array to be filled and fills it according to rows as specified in array of colums
*/
{


  int errcode;
  

    if(debug)
    {
      for(int i=0; i<rowsize; i++)
	{
	  bootrows[i]=i;
	}
    }
    

	
    //Let's fill in bootstnp
    for(int i=0; i<rowsize; i++)
      {
	errcode = bootrows[i];
	
	cblas_dcopy(colsize,&omatrix[index(errcode,0,colsize)],1,&bootmatrix[index(i,0,colsize)],1);
      }
      
   
}

void printProgBar(int percent){
  char bar[50];
  
  for(int i=0; i<50; i++)
    {
      if(i<(percent/2))
	{
	  bar[i]='=';

	}
      else 
	{
	  if (i==(percent/2))
	    {
	      bar[i]='>';

	    }
	  else
	    {
	      bar[i]=' ';
	    }
	}
    }
  cout<<"\r" "["<<bar<<"] ";
  cout.width(3);
  cout<<percent<<"%     "<<flush;
}

void MatSlice(double *&slicemat, double *&omat,int colstart,int colsize,int rowsize, int zsize)
//return 2d slices of a 3d array
{
  for(int i=0; i<zsize; i++)
    {

      cblas_dcopy(rowsize,&omat[multindex(0,colstart,i,rowsize,colsize,zsize)],colsize,&slicemat[index(0,i,zsize)],zsize);
    }
}


herr_t GetSlice(int rowstart,int colstart,int bsistart,int rownum, int colnum,int bsinum,double *&matrix,hid_t file, hid_t dataset,hid_t dataspace)
{

  //This in theory is the dimension of the ouput data (All snps, one gene, all bootstrap iterations)
  hsize_t slabsize[3]={rownum,colnum,bsinum};

  hsize_t offset[3]={rowstart,colstart,bsistart};

  hsize_t offset_out[3]={0,0,0};
  hid_t memspace;
  herr_t status;


  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,slabsize,NULL);
  
  memspace = H5Screate_simple(3,slabsize,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,slabsize,NULL);
  status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,matrix);

   H5Sclose(memspace);
   return(status);
}


herr_t ReadMatrix(double *&matrix,int rtot,int ctot,int rstart,int rsize, int cstart, int csize,const char* filename,const char* snpgene,const char*lockfile)
{
  struct stat buffer;
  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;
  FILE * fp;
  hsize_t dimsout[2]={rsize,csize};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={rsize,csize};

   
  
  matrix = (double *)mkl_malloc(rsize*csize*sizeof(double),64);
  
  //  cout<<"obtaining file lock"<<endl;
  while(stat(lockfile,&buffer)==0)
    {
      usleep(999999);
      cout<<lockfile<<" exists!"<<endl;
    }
  fp = fopen(lockfile,"ab+");
  
  //  cout<<"file lock obtained"<<endl;

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

  fclose(fp);
  remove(lockfile);

  return(status);




}

void ReadAnno (char **&name, int *&array, int numstart,int numsize,const char* filename,int arraydim)
{
  
  hid_t dataspacename,dataspacearray,file_id,datasetname,datasetarray,rank,memspacename,memspacearray;
  hid_t fstrtype;
  hid_t mstrtype;
  hsize_t offsetname[1]={numstart};
  hsize_t offsetarray[2]={numstart,0};
  hsize_t dimsarray[2]={numsize,arraydim};
  hsize_t dimsname[1]={numsize};
  hsize_t offset_outname[1]={0};
  hsize_t offset_outarray[2]={0,0};
  herr_t status;
  
  
  array = (int *)mkl_malloc(arraydim*numsize*sizeof(int),64);
  name = (char**)mkl_malloc(numsize*sizeof(char*),64);
  
  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(arraydim==2)
    {
      datasetname = H5Dopen2(file_id,"snpname",H5P_DEFAULT);
      datasetarray = H5Dopen2(file_id,"snparray",H5P_DEFAULT);
      
    }
  else
    {
      datasetname = H5Dopen2(file_id,"genename",H5P_DEFAULT);
      datasetarray = H5Dopen2(file_id,"genearray",H5P_DEFAULT);
    }
  dataspacename = H5Dget_space(datasetname);
  dataspacearray = H5Dget_space(datasetarray);
  
  status = H5Sselect_hyperslab(dataspacename,H5S_SELECT_SET,offsetname,NULL,dimsname,NULL);
  status = H5Sselect_hyperslab(dataspacearray,H5S_SELECT_SET,offsetarray,NULL,dimsarray,NULL);

  memspacename=H5Screate_simple(1,dimsname,NULL);
  memspacearray=H5Screate_simple(2,dimsarray,NULL);
  

  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  mstrtype = H5Tcopy(H5T_C_S1);
  status= H5Tset_size(mstrtype,H5T_VARIABLE);
  
  status = H5Sselect_hyperslab(memspacename,H5S_SELECT_SET,offset_outname,NULL,dimsname,NULL);
  status= H5Sselect_hyperslab(memspacearray,H5S_SELECT_SET,offset_outarray,NULL,dimsarray,NULL);

  status = H5Dread(datasetname,fstrtype,memspacename,dataspacename,H5P_DEFAULT,name);
  status = H5Dread(datasetarray,H5T_NATIVE_INT,memspacearray,dataspacearray,H5P_DEFAULT,array);

  status = H5Tclose(fstrtype);
  status = H5Tclose(mstrtype);
  status = H5Dclose(datasetname);
  status = H5Dclose(datasetarray);
  status = H5Sclose(dataspacename);
  status = H5Sclose(dataspacearray);
  status = H5Fclose(file_id);
  

}



void GetQuantile(double *&matrix, size_t rowsize,size_t colsize, size_t bootsize,double *&quantiles)
{

  int deciles = 2;
  VSLSSTaskPtr task;
  int status;
  double variation[colsize*rowsize];
  double mean[colsize*rowsize];
  double r2[colsize*rowsize];
  double o_quant[2]={0.025,0.975};

  double o_stat[colsize*rowsize*bootsize];
  MKL_INT x_storage =VSL_SS_MATRIX_STORAGE_COLS;
  MKL_INT o_storage = VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT64 estimate=0;
  int p = rowsize*colsize;
  int n=bootsize;
  

  status = vsldSSNewTask( &task,&p,&n,&x_storage,(double*)matrix,NULL,NULL);


  status = vsldSSEditQuantiles(task,&deciles,o_quant,(double*)quantiles,(double*)o_stat,&o_storage);

  status = vsldSSCompute(task,VSL_SS_QUANTS,VSL_SS_METHOD_FAST);

  status = vslSSDeleteTask(&task);




  
}

double* TestGenerate(double *&matrix, int rowsize,int colsize,int rowstart, int colstart,int zindex)
//Create a 2d test array for the rest of the functions
{ 
  for(int i=0; i<rowsize;i++)
    {
      for(int j=0; j<colsize;j++)
	{
	  matrix[index(i,j,colsize)]=((double) rowstart+(double)i)+((double) colstart+(double)j*0.1)+(double)zindex*0.01;
	}
    }
  return(matrix);
}

bool isCis(const int* snparray,const int* genearray,int snp,int gene,int cisdist)
{
  int snpchr = snparray[index(snp,0,2)];
  int snppos = snparray[index(snp,1,2)];
  int genechr = genearray[index(gene,0,3)];
  int genestart = genearray[index(gene,1,3)];
  int geneend = genearray[index(gene,2,3)];
  
  return( (snpchr==genechr)&&((abs(snppos-genestart)<cisdist)||(abs(snppos-geneend)<cisdist)));
  
}
void CisTransOut (const double *quantilearray,const int snpstart,const int snpsize,const int genestart,const int genesize, const int casesize,const int zsize,const char *annofile,double t_thresh,const char*cisfile,const char* transfile,int cisdist,const char* lockfile,int mychunk,const char* progfile)
  
{

  struct stat buffer;
  int status;
  FILE* fp;
  ofstream cis,trans,prog;
  double cor1,cor2,worstcor;
  double cordenom;
  int *snparray,*genearray;
  char** snpnames,**genenames;

  while(stat(lockfile,&buffer)==0)
    {
      usleep(9999);
      cout<<lockfile<<" exists!"<<endl;
    }
  fp = fopen(lockfile,"ab+");

  ReadAnno(snpnames,snparray,snpstart,snpsize,annofile,2);
  ReadAnno(genenames,genearray,genestart,genesize,annofile,3);
  
  cordenom = sqrt(casesize-2);


  prog.open(progfile,ofstream::out|ofstream::app);
  prog<<mychunk<<endl;
  cis.open(cisfile,ofstream::out|ofstream::app);
  trans.open(transfile,ofstream::out|ofstream::app);
  if(mychunk==0){
    cis<<"SNP\tGene\tt-stat1\tt-stat2"<<endl;
    trans<<"SNP\tGene\tt-stat1\tt-stat2"<<endl;
  }
  
  for(int i=0; i<snpsize; i++)
    {
      for(int j=0; j<genesize; j++)
	{
	  cor1=quantilearray[multindex(i,j,0,snpsize,genesize,zsize)];
	  cor1=cordenom*(cor1/sqrt(1-(cor1*cor1)));
	  cor2=quantilearray[multindex(i,j,1,snpsize,genesize,zsize)];
	  cor2=cordenom*(cor2/sqrt(1-(cor2*cor2)));
	  worstcor = fabs(cor1)>fabs(cor2) ? cor2 : cor1;
	  //	  cout<<snpnames[i]<<"-"<<genenames[j]<<": "<<fabs(worstcor)<<" "<<t_thresh<<endl;
	  

	  if(fabs(worstcor)>t_thresh && ((cor1<0)==(cor2<0))&&((!isnan(cor1))&&(!isnan(cor2))))
	    {
	      if(isCis(snparray,genearray,i,j,cisdist))
		{
		  cis<<snpnames[i]<<"\t"<<genenames[j]<<"\t"<<cor1<<"\t"<<cor2<<endl;
		}
	      else
		{
		  trans<<snpnames[i]<<"\t"<<genenames[j]<<"\t"<<cor1<<"\t"<<cor2<<endl;
		}
	    }
	}
    } 

  prog.close();
  trans.close();
  cis.close();
  fclose(fp);
  remove(lockfile);

}






void readparam(hid_t file_id, const char* attributename,char*&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;

  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  fstrtype= H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,fstrtype,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}
	  
void readparam(hid_t file_id,const char* attributename,int&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,H5T_NATIVE_INT,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}

void readparam(hid_t file_id,const char* attributename,double&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,H5T_NATIVE_DOUBLE,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}
#endif
