#ifndef MKTEST_H
#define MKTEST_H
#include "hdf5.h"
#include "hdf5_hl.h"

enum snp_exp {SNP,EXP};

int index (int i, int j, int M)
{
  return((M*i)+j);
}

void scale(double matrix[],int r,int c)
{

  double colmean=0.0;
  double colsd=0.0;
  for(int j=0;j<c;j++)
    {
      colmean=0;
      colsd=0;
      for(int i=0;i<r;i++)
	{
	  colmean+=matrix[index(i,j,c)];
	}
      colmean = colmean/r;

      for(int i=0;i<r;i++)
	{
	  colsd+=pow(matrix[index(i,j,c)]-colmean,2);
	}
      colsd=sqrt(colsd/(r-1));
      for(int i=0;i<r;i++)
	{
	  if(colsd==0)
	    {
	      matrix[index(i,j,c)]=0;
	    }
	  else
	    {
	      matrix[index(i,j,c)]=(matrix[index(i,j,c)]-colmean)/(colsd);
	    }
	}
      //	cout<<"\n";
    }

}

void bootstrap(double *bootsnpmatrix,double *bootexpmatrix, const double osnpmatrix[], const double oexpmatrix[], int r, int m, int n,int bs)
{
  int *bootcols;
  

  bootcols = (int *)mkl_malloc(r*sizeof(int),64);
  if(bs==0){
    bootsnpmatrix= (double *) mkl_malloc(r*m*sizeof(double),64);
    bootexpmatrix= (double *) mkl_malloc(r*n*sizeof(double),64);
  }
  
  for(int i=0; i<r; i++)
      {
	bootcols[i]=rand()%r;
	
      }
    
    //Let's fill in bootstnp
    for(int i=0; i<r; i++)
      {
	for(int j=0; j<m; j++)
	  {
	    bootsnpmatrix[index(i,j,m)] = osnpmatrix[index(bootcols[i],j,m)];
	  }
	for(int k=0; k<n; k++)
	  {
	    bootexpmatrix[index(i,k,n)] = oexpmatrix[index(bootcols[i],k,n)];
	  }
      }
    mkl_free(bootcols);
}

void readmatrix(double *matrix,int rtot,int ctot,int r,int c, const char* filename,snp_exp snpgene)
{
  hid_t dataspace,file,dataset,dataspace,rank,status_n,memspace;
  herr_t status;
  hsize_t dimsout[2]={rtot,ctot};
  hsize_t offset[2]={0,0};
  hsize_t dimsm[2] ={r,c};
  
  matrix = (double *)mkl_malloc(r*c*sizeof(double),64);
  file  = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if(snpgene==SNP)
    {
      dataset = H5Dopen(file,"snps",H5P_DEFAULT);
    }else
    {
      dataset = H5Dopen(file,"genes",H5P_DEFAULT);
    }
  dataspace = H5Dget_space(dataset);
  rank = H5Sget_simple_extent_ndims(dataspace);
  status_n = H5Sget_simple_extent_dims(dataspace,dimsout,NULL);
  


}

void readanno (string name[], int chr[], int posstart[], int end[], int num, const char* filename,snp_exp snpgene)
{

  name = new string[num];
  chr = (int *)mkl_malloc(num*sizeof(int),64);
  posstart = (int *)mkl_malloc(num*sizeof(int),64);
  if(snpgene!=SNP){
    end = (int *) mkl_malloc(num*sizeof(int),64);
  }



  ifstream annofile;
  annofile.open(filename);
  for(int i=0;i<num;i++)
    {
      if(snpgene==SNP){
	annofile>>name[i]>>chr[i]>>posstart[i];
      }else
	annofile>>name[i]>>chr[i]>>posstart[i]>>end[i];
    }
  annofile.close();
}



void writemat(double matrix[],int r, int c, hid_t filespace,hid_t dset_id,int bsindex)
{
  hid_t plist_xfer;
  hsize_t slabsize[3]={r,c,1};
  hsize_t offset[3]={0,0,bsindex};
  hid_t memspace = H5Screate_simple(3,slabsize,NULL);
  herr_t status;

  filespace = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,slabsize,NULL);

  plist_xfer=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_buffer(plist_xfer,(hsize_t)r*(hsize_t)c,NULL,NULL);
  status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_xfer,matrix);
}	  

#endif
