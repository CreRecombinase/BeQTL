#ifndef MKTEST_HPP
#define MKTEST_HPP
#include <math.h>
#include "mkl.h"
#include "strideiter.hpp"

#include <algorithm>



bool debug=false;
using namespace std;



void printProgBar(int percent){
  cout<<"[";
  for(int i=0; i<50; i++)
    {
      if(i<(percent/2))
	{
	  cout<<"=";

	}
      else 
	{
	  if (i==(percent/2))
	    {
	      cout<<">";

	    }
	  else
	    {
	      cout<<" ";
	    }
	}
    }
  cout<<"]"<<percent<<" %\r";
  cout.flush();
}

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




void GetQuantile(double *&matrix, size_t rowsize,size_t colsize, size_t bootsize,double *&quantiles)
//Let's try to do the quantile estimation in four chunks if we can't get quantile estimation to work off the bat
{

  int qsize=1;
  int deciles = 3;
  VSLSSTaskPtr task;
  int status;
  double variation[colsize*rowsize];
  double mean[colsize*rowsize];
  double r2[colsize*rowsize];
  double o_quant[3]={0.025,0.5,0.975};

  double* o_stat=NULL;
  double*tmatrix = NULL;
  MKL_INT x_storage =VSL_SS_MATRIX_STORAGE_COLS;
  MKL_INT o_storage = VSL_SS_MATRIX_STORAGE_ROWS;
  unsigned MKL_INT64 estimate=0;
  int p = rowsize*colsize;
  int n=bootsize;

  o_stat = (double*)mkl_malloc(colsize*rowsize*bootsize*sizeof(double),64);
  if(o_stat==NULL||debug)
    {
      qsize=4;
      p=(rowsize/qsize)*colsize;
      if(o_stat==NULL){
	o_stat=(double*)mkl_malloc(colsize*(rowsize/qsize)*bootsize*sizeof(double),64);
      }
      tmatrix=(double*)mkl_malloc(colsize*(rowsize/qsize)*bootsize*sizeof(double),64);
      for(int q=0; q<qsize; q++)
	{
	  cout<<"q iteration "<<q<<endl;
	  for(int i=0; i<(rowsize/qsize);i++)
	    {
	      for(int j=0;j<colsize;j++)
		{
		  cblas_dcopy(bootsize,&matrix[index((q*(rowsize/qsize))+i,j,colsize)],colsize*rowsize,&tmatrix[index(i,j,colsize)],colsize*(rowsize/qsize));
		}
	    }

	  
	  status = vsldSSNewTask( &task,&p,&n,&x_storage,(double*)tmatrix,NULL,NULL);
	  status = vsldSSEditQuantiles(task,&deciles,o_quant,(double*)&quantiles[multindex(q*(rowsize/qsize),0,0,rowsize,colsize,2)],(double*)o_stat,&o_storage);
	  status = vsldSSCompute(task,VSL_SS_QUANTS,VSL_SS_METHOD_FAST);
	  status = vslSSDeleteTask(&task);
	  
	}
      mkl_free(o_stat);
      mkl_free(tmatrix);
    }
  else{

    status = vsldSSNewTask( &task,&p,&n,&x_storage,(double*)matrix,NULL,NULL);


    status = vsldSSEditQuantiles(task,&deciles,o_quant,(double*)quantiles,(double*)o_stat,&o_storage);

    status = vsldSSCompute(task,VSL_SS_QUANTS,VSL_SS_METHOD_FAST);

    status = vslSSDeleteTask(&task);
    mkl_free(o_stat);
  }



  
}

void TestGenerate(double *matrix, size_t rowsize,size_t colsize,size_t rowstart, size_t colstart,size_t zindex)
//Create a 2d test array for the rest of the functions
{ 
  for(int i=0; i<rowsize;i++)
    {
      for(int j=0; j<colsize;j++)
	{
	  matrix[index(i,j,colsize)]=((double) rowstart+(double)i)+((double) colstart+(double)j*0.1)+(double)zindex*0.01;
	}
    }

}
#endif
