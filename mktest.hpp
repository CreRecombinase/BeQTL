#ifndef MKTEST_HPP
#define MKTEST_HPP
#include <string>
#include "hdf5.h"
#include <math.h>
#include "mkl.h"
#include "hdf5_hl.h"
#include <sys/file.h>
#include <omp.h>

bool debug=false;
using namespace std;
enum snp_exp {SNP,GENE};

int multindex(int row,int column,int zindex,int nrows,int ncols)
{
  return((nrows*ncols*zindex)+(ncols*row)+column);
}

int index (int row, int column, int nocols)
//Function for indexing 2 dimensional, row major order array
//row     row index
//column  column index
//nocols  number of columns
{
 
  return((nocols*row)+column);
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
  if(debug)
    {
      PrintMat(means,1,cols);
      PrintMat(variances,1,cols);
    }
  for(int i=0; i<rows; i++)
    {
  vdSub(cols,&matrix[index(i,0,cols)],means,&matrix[index(i,0,cols)]);
  vdDiv(cols,&matrix[index(i,0,cols)],variances,&matrix[index(i,0,cols)]);
    }

}

void MakeBootRows(int *&bootrows, int rowsize,int bstotal)
//Function to geneerate 2d array of bootstrap rows
{
  VSLStreamStatePtr stream;
  int errcode;
  
  bootrows = (int*) mkl_malloc(rowsize*bstotal*sizeof(int),64);
  errcode = vslNewStream(&stream,VSL_BRNG_MCG31,123);
  errcode = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,rowsize*bstotal,bootrows,0,rowsize);

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
	cblas_dcopy(colsize,&omatrix[index(bootrows[i],0,colsize)],1,&bootmatrix[index(i,0,colsize)],1);
      }
   
}

void MatSlice(double *&slicemat, double *&omat,int colstart,int colsize,int rowsize, int zsize)
//return 2d slices of a 3d array
{
  for(int i=0; i<zsize; i++)
    {
      //      cout<<omp_get_thread_num()<<endl;
      cblas_dcopy(rowsize,&omat[multindex(0,colstart,i,rowsize,colsize)],colsize,&slicemat[index(0,i,zsize)],zsize);
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

  //  matrix = (double *) mkl_malloc(rownum*colnum*bsinum*sizeof(double),64);


  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,slabsize,NULL);
  
  memspace = H5Screate_simple(3,slabsize,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,slabsize,NULL);
  status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,matrix);

   H5Sclose(memspace);
   return(status);
}


herr_t ReadMatrix(double *&matrix,int rtot,int ctot,int rstart,int rsize, int cstart, int csize,const char* filename,const char* snpgene,const char*lockfile)
{
  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;
  int result=-1;
  int fd=-1;
  hsize_t dimsout[2]={rsize,csize};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={rsize,csize};
  
  matrix = (double *)mkl_malloc(rsize*csize*sizeof(double),64);
  while(fd<0|result<0)
    {
      fd = open(lockfile,O_RDONLY);

      result = flock(fd,LOCK_EX);
    }

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
  
  result = flock(fd,LOCK_UN);
  close(fd);
  return(status);




}

void ReadAnno (string name[], int array[], int numstart,int numend, const char* filename,int arraydim)
{
  
  int numsize = (numend-numstart)+1;
  hid_t dataspacename,dataspacearray,file_id,datasetname,datasetarray,rank,memspacename,memspacearray;
  hid_t fstrtype;
  hid_t mstrtype;
  name = new string[numsize];
  hsize_t offsetname[1]={numstart};
  hsize_t offsetarray[2]={numstart,0};
  hsize_t dimsarray[2]={numsize,arraydim};
  hsize_t dimsname[1]={numsize};
  hsize_t offset_outname[1]={0};
  hsize_t offset_outarray[2]={0,0};
  herr_t status;
  
  
  array = (int *)malloc(arraydim*numsize*sizeof(int));
  
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

herr_t WriteMat(double matrix[],int rowstart,int rowsize, int colstart,int colsize,int rowtotal,int coltotal,const char* filename,int snpchunk,const char* writefilelock)
{

  hid_t file,filespace,dataset;
  hsize_t space_dims[3]={rowtotal,coltotal,2};
  int wlockf=-1;
  int result=-1;
  hid_t plist_xfer;
  hsize_t slabsize[3]={rowsize,colsize,2};
  hsize_t offset[3]={rowstart,colstart,0};
  hid_t memspace = H5Screate_simple(3,slabsize,NULL);
  herr_t status;

  while(wlockf<0|result<0)
    {
      wlockf=open(writefilelock,O_RDONLY);
      result=flock(wlockf,LOCK_EX);
    }
  
  if(snpchunk==0)
    {
      file=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      filespace = H5Screate_simple(3,space_dims,NULL); 
      dataset= H5Dcreate(file,"quantiles",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    }
  else{
    file = H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
    dataset = H5Dopen(file,"quantiles",H5P_DEFAULT);
    filespace = H5Dget_space(dataset);
  }


  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,slabsize,NULL);


  plist_xfer=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist_xfer,(hsize_t)rowsize*(hsize_t)colsize,NULL,NULL);
  status=H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace,filespace,plist_xfer,matrix);
  H5Sclose(memspace);

  status = H5Dclose(dataset);
  status = H5Sclose(filespace);
  status = H5Fclose(file);  

  result=flock(wlockf,LOCK_UN);
  result=close(wlockf);
  return(status);

}



void GetQuantile(double *&matrix, int rowsize,int colsize, int bootsize,double *&quantiles)
{

  int deciles = 2;
  VSLSSTaskPtr task;
  int status;
  double variation[colsize*rowsize];
  double mean[colsize*rowsize];
  double r2[colsize*rowsize];
  double o_quant[2]={0.05,0.95};
  //  double quantiles[cols*2];
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

void TestGenerate(double *matrix, int rowsize,int colsize,int colstart,int zindex)
//Create a 2d test array for the rest of the functions
{ 
  for(int i=0; i<rowsize;i++)
    {
      for(int j=0; j<colsize;j++)
	{
	  matrix[index(i,j,colsize)]=(double)i+((double) colstart+(double)j*0.1)+(double)zindex*0.01;
	}
    }
}

bool isCis(const double* snparray,const double* genearray,int snp,int gene,int cisdist)
{
  int snpchr = snparray[index(snp,0,2)];
  int snppos = snparray[index(snp,1,2)];
  int genechr = genearray[index(gene,0,3)];
  int genestart = genearray[index(gene,1,3)];
  int geneend = genearray[index(gene,2,3)];
  
  return( (snpchr==genechr)&((abs(snppos-genestart)<cisdist)|(abs(snppos-geneend)<cisdist)));
  
}
void cistransout (const double *tmatrix,int snpsize,int genesize, const double* snparray,const double* genearray, const string  *snpnames,const string *genenames,double t_thresh,const char*cisfile,const char* transfile,int cisdist)
  
{
  FILE* cis,*trans;
  
  cis = fopen(cisfile,"w");
  trans = fopen(transfile,"w");
  
  for(int i=0; i<snpsize; i++)
    {
      for(int j=0; j<genesize; j++)
	{
	  if(abs(tmatrix[index(i,j,genesize)])>t_thresh)
	    {
	      if(isCis(snparray,genearray,i,j,cisdist))
		{
		  cis<<snpnames[i]<<"\t"<<genenames[j]<<"\t"<<tmatrix[index(i,j,genesize)];
		}
	      else
		{
		  trans<<snpnames[i]<<"\t"<<genenames[j]<<"\t"<<tmatrix[index(i,j,genesize)];
		}
	    }
	}
    }	  

  fclose(trans);
  fclose(cis);
}


void worstT (double *&matrix,const double* corarray,int rowsize,int colsize,int samplesize,double *tmat)		
//Calcuate t stat for correlation bootstrap and return the 'worst' of the two t-stats
{

  double tempt1,tempt2;
  double tcorr;
  double tr;
  tcorr = sqrt(samplesize-2);
  matrix = (double*)mkl_malloc(rowsize*colsize*sizeof(double),64);
  
  for(int i=0; i<rowsize; i++)
    {
      for(int j=0; j<colsize; j++)
	{
	  tr = corarray[multindex(i,j,0,rowsize,colsize)];
	  tempt1=tcorr*(tr/sqrt(1-(tr*tr)));
	  tr = corarray[multindex(i,j,1,rowsize,colsize)];
	  tempt2=tcorr*(tr/sqrt(1-(tr*tr)));
	  if(tempt1>0&tempt2>0)
	    {
	      if(tempt1>tempt2)
		{
		  matrix[index(i,j,colsize)]=tempt2;
		}
	      else
		{
		  matrix[index(i,j,colsize)]=tempt1;
		}
	    }
	  else
	    {
	      if(tempt1<0&tempt2<0)
		{
		  if(tempt1>tempt2)
		    {
		      matrix[index(i,j,colsize)]=tempt1;
		    }
		  else
		    {
		      matrix[index(i,j,colsize)]=tempt2;
		    }
		}
	      else
		{
		  matrix[index(i,j,colsize)]=0;
		}
	    }
	}
    }
}

#endif
