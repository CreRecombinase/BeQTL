#include "mkl.h"
#include <stdio.h>
#include <iostream>
using namespace std;
#define DIM 3  //Dimension of the task
#define N 4 // Number of observations
#define M 100 //accuracy of quantile computation
#define EPS 0.01

int main()
{
  int i, status,j;
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
  q_order_n = 3;
  xstorage = VSL_SS_MATRIX_STORAGE_ROWS;
  params = EPS;
  nparams = VSL_SS_SQUANTS_ZW_PARAMS_N;
  for(i=0; i<DIM; i++){
    for(j=0; j<N; j++){
      x[i*N+j]= i*N+j;
      cout<<x[i*N+j]<<"\t";
    }
    cout<<"\n";
  }
  cout<<"---------\nq_order:"<<endl;
  // Calculate percentiles
  q_order[0] = 0.025;
  q_order[1] = 0.5;
  q_order[2] = 0.975;
  for( i=0; i< M; i++){
    cout<<q_order[i]<<"\t";
  }
  cout<<"\n-----"<<endl;
  // Create task
  status = vslsSSNewTask( &task, &p, &n, &xstorage, x, 0, indices);

  status = vslsSSEditStreamQuantiles(task,&q_order_n,q_order,
				     quants,&nparams,&params);

  //Compute percentile with accuracy eps

  status = vslsSSCompute (task,VSL_SS_STREAM_QUANTS, 
			  VSL_SS_METHOD_SQUANTS_ZW);

  //Deallocate the task resources 
  status = vslSSDeleteTask( &task);
  for(i=0; i<M; i++){
    cout<<quants[i]<<"\t";
  }
  
  return(0);
}



