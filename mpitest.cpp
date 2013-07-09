#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include "mkl_blacs.h"

#include "mpi.h"

using namespace std;
  
int index (int i, int j, int M)
{
  return((M*i)+j);
}

int main(int argc,char **argv)
{
  int i,j,k;
  // MPI
  int myrank_mpi, nprocs_mpi;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs_mpi);
  //BLACS
  
  int ictxt,nprow,npcol,myrow_snp,myrow_exp,myrow_cor,mycol_snp,mycol_exp,mycol_cor,nb;
  int ZERO=0,ONE=1;
  int info,itemp;

  nprow=3;npcol=3;nb=2;

  blacs_pinfo_(&myrank_mpi,&nprocs_mpi);
  blacs_get_(-1,0,&ictxt);
  blacs_gridinit_(&ictxt,"Row",nprow,npcol);
  blacs_gridinfo_(ictxt,&nprow,&npcol,&myrow,&mycol);

  int rows = 20;
  int snps = 30;
  int genes = 25;
  double *snpmat;
  double *expmat;
  double *cormat;

  int procrows=3;
  int proccols=3;
  
  int snpdimensions[4];
  int expdimensions[4];
  int cordimensions[4];

  //rows of global matrix
  snpdimensions[0]=rows;
  expdimensions[0]=rows;
  cordimensions[0]=snps;
  //columns of global matrix
  snpdimensions[1]=snps;
  expdimensions[1]=genes;
  cordimensions[1]=genes;
  //rows of local matrix
  snpdimensions[2]=4;
  expdimensions[2]=4;
  cordimensions[2]=4;
  //columns of local matrix
  snpdimensions[3]=2;
  expdimensions[3]=2;
  cordimensions[3]=2;

  //compute how many rows each process gets
  int nrows_snp = numroc_(&rows,&snpdimensions[2],&myrow_snp,&iZERO,&procrows);
  int nrows_gene = numroc_(&rows,&expdimensions[2],&myrow_exp,&iZERO,&procrows);
  int nrows_cor = numroc(&snps, &cordimensions[2],&myrow_cor,&iZERO,&procrows);

  int ncols_snp = numroc_(&snps,&snpdimensions[3],&mycol_snp,&iZERO,&proccols);
  int ncols_exp = numroc_(&genes,&expdimensions[3],&mycol_exp,&iZERO,&proccols);
  int ncols_cor = numroc_(&genes,&cordimensions[3],&mycol_cor,&iZERO,&proccols);


  //fill matrix description arrays for BLACS matrices
  int descsnp[9],descexp[9],desccor[9];
  //descinit_(descsnp,&rows,&snps,&nb,&nb,&iZERO,&iZERO,&ictxt,&nrows_snp,&info);
  //descinit_(descexo,&rows,&genes,&nb,&nb,&iZERO,&iZERO,&ictxt,&nrows_exp,&info);
  //descinit_(desccor,&snps,&genes,&nb,&nb,&iZERO,&iZERO,&ictxt,&nrows_cor,&info);

  //declare local memory for submatrices
  //snpmat = (double*) malloc(nrows_snp*ncols_snp*(sizeof(double)));
  //expmat = (double*) malloc(nrows_exp*ncols_exp*(sizeof(double)));
  //snpmat = (double*) malloc(nrows_cor*ncols_cor*(sizeof(double)));

  //Fill submatrices here using hdf5

  cout<<"Rank:"<<myrank_mpi<<"\tnrows_snp:"<<nrows_snp<<"\tncols_snp:"<<ncols_snp<<endl;

  Cblacs_gridexit(0);
  MPI_Finalize();
  return(0);
}

  
  
 
    
  

  
