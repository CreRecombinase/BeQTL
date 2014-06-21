#include "mkl.h"
#include <stdio.h>

#define DIM 3  //Dimension of the task
#define N 1000 // Number of observations
#define M 100 //accuracy of quantile computation
#define EPS 0.01

int main()
{
  int i, status;
  VSLSSTaskPtr task;
  float x[DIM*N]; //matrix of observations
  float q_order[M],quants[M];
  float params;
  MKL_INT q_order_n;
  MKL_INT p,n,nparams,xstorage;
  int indices[DIM]={1,0,0}; //1st vector component will be processed
  //Parameters of the task and initialization

  p = DIM;
  n = N;
  q_order_n = M;
  xstorage = VSL_SS_MATRIX_STORAGE_ROWS;
  params = EPS;
  nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
 
  // Calculate percentiles
  for( i=0; i< M; i++){
    q_order[i] = (float)i/(float)M;
  }
  
  // Create task
  status = vslsSSNewTask( &task, &p, &n, &xstorage, x, 0, indices);

  status = vslsSSEditStreamQuantiles(task,&q_order_n,q_order,
				     quants,&nparams,&params);

  //Compute percentile with accuracy eps

  status = vslsSSCompute (task,VSL_SS_STREAM_QUANTS, 
			  VSL_SS_METHOD_SQUANTS_ZW);

  //Deallocate the task resources 
  status = vslSSDeleteTask( &task);
  
  return(0);
}



