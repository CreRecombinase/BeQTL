#ifndef MKTEST_HPP
#define MKTEST_HPP
#include <string>
#include "hdf5.h"
#include <math.h>
#include "mkl.h"
#include "hdf5_hl.h"

bool debug=false;
using namespace std;
enum snp_exp {SNP,GENE};

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

void GetSlice(const char* filename,int rowstart,int colstart,int bsistart,int rownum, int colnum,int bsinum,double *&matrix)
{
  hid_t file;
  hsize_t slabsize[3]={rownum,colnum,bsinum};
  hsize_t offset[3]={rowstart,colstart,bsistart};
  //This in theory is the dimension of the ouput data (All snps, one gene, all bootstrap iterations)
  hsize_t dimsm[3]={rownum,colnum,bsinum};
  hsize_t offset_out[3]={0,0,0};
  hid_t dataset,dataspace,filespace,memspace;
  herr_t status;

  //  matrix = (double *) mkl_malloc(rownum*colnum*bsinum*sizeof(double),64);
  file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  dataset = H5Dopen2(file,"cor",H5P_DEFAULT);
  dataspace = H5Dget_space(dataset);

  status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);
  
  memspace = H5Screate_simple(3,dimsm,NULL);

  status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);
  status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,matrix);
  status= H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);
}


void ReadMatrix(double *&matrix,int rtot,int ctot,int rstart,int rend, int cstart, int cend,const char* filename,const char* snpgene)
{
  hid_t dataspace,file,dataset,rank,memspace;
  herr_t status;
  int rsize = (rend-rstart)+1;
  int csize = (cend-cstart)+1;
  hsize_t dimsout[2]={rsize,csize};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={rsize,csize};
  
  matrix = (double *)mkl_malloc(rsize*csize*sizeof(double),64);
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


void InitDataFile(const char* filename,const int rowtotal,const int coltotal,const int ztotal,const char* dataname)
{
  hid_t file_id,dset_id,filespace;
  hsize_t cordims[3];

  file_id = H5Fcreate (filename,H5F_ACC_RDWR,H5P_DEFAULT,H5P_DEFAULT);



  cordims[0]=rowtotal;
  cordims[1]=coltotal;
  cordims[2]=ztotal;

  filespace = H5Screate_simple(3,cordims,NULL);
  
  dset_id = H5Dcreate2(file_id,dataname,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Sclose(filespace);
}
void WriteMat(double matrix[],int rowstart,int rowend, int colstart,int colend,int bsindex,const hid_t dataset, const hid_t filespace)
{

  int rowsize = (rowend-rowstart)+1;
  int colsize = (colend-colstart)+1;
  
  hid_t plist_xfer;
  hsize_t slabsize[3]={rowsize,colsize,1};
  hsize_t offset[3]={rowstart,colstart,bsindex};
  hid_t memspace = H5Screate_simple(3,slabsize,NULL);
  herr_t status;

  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,slabsize,NULL);


  plist_xfer=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist_xfer,(hsize_t)rowsize*(hsize_t)colsize,NULL,NULL);
  status=H5Dwrite(dataset,H5T_NATIVE_DOUBLE,memspace,filespace,plist_xfer,matrix);
  status = H5Sclose(memspace);
}



 void GetQuantile(double *&matrix, int rows, int cols,double *&quantiles)
{

  int deciles = 2;
  VSLSSTaskPtr task;
  int status;
  double variation[cols];
  double mean[cols];
  double min_est[cols];
  double max_est[cols];
  double r2[cols];
  double o_quant[2]={0.05,0.95};
  //  double quantiles[cols*2];
  double o_stat[cols*rows];
  MKL_INT x_storage =VSL_SS_MATRIX_STORAGE_ROWS;
  MKL_INT o_storage = VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT64 estimate=0;
  int dim=1;
  int n=rows*cols;

  status = vsldSSNewTask( &task,&cols,&rows,&x_storage,(double*)matrix,NULL,NULL);
 

  status = vsldSSEditQuantiles(task,&deciles,o_quant,(double*)quantiles,(double*)o_stat,&o_storage);
  status = vsldSSCompute(task,VSL_SS_QUANTS,VSL_SS_METHOD_FAST);


  
}

void TestGenerate(double *&matrix, int rowsize,int colsize,int zindex)
//Create a 2d test array for the rest of the functions
{
  matrix = (double*)mkl_malloc(rowsize*colsize*sizeof(double),64);
  
  for(int i=0; i<rowsize;i++)
    {
      for(int j=0; j<colsize;j++)
	{
	  matrix[index(i,j,colsize)]=(double)i+(double)j*0.1+(double)zindex*0.01;
	}
    }
}

#endif
