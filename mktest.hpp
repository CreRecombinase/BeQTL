#ifndef MKTEST_H
#define MKTEST_H
#include "hdf5.h"
#include "hdf5_hl.h"

enum snp_exp {SNP,GENE};

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

void readmatrix(double *matrix,int rtot,int ctot,int rstart,int rend, int cstart, int cend,const char* filename,snp_exp snpgene)
{
  hid_t dataspace,file,dataset,dataspace,rank,status_n,memspace;
  herr_t status;
  hsize_t dimsout[2]={(rend-rstart)+1,(cend-cstart)+1};
  hsize_t offset[2]={rstart,cstart};
  hsize_t offset_out[2] ={0,0};
  hsize_t dimsm[2] ={(rend-rstart)+1,(cend-cstart)+1};
  
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

  status = H5select_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,dimsm,NULL);

  memspace = H5Screate_simple(RANK_OUT,dimsm,NULL);
  status = H5select_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,dimsm,NULL);
  status = H5Dread(dataset,H5T_NATIVE_DOUBLE,memspace,dataspace,H5P_DEFAULT,matrix);
  


}

void readanno (string name[], int chr[], int posstart[], int end[], int numstart,int numend, const char* filename,snp_exp snpgene)
{

  int numsize = (numend-numstart)+1
  name = new string[numsize];
  chr = (int *)mkl_malloc(numsize*sizeof(int),64);
  posstart = (int *)mkl_malloc(numsize*sizeof(int),64);
  if(snpgene!=SNP){
    end = (int *) mkl_malloc(numsize*sizeof(int),64);
  }



  ifstream annofile;
  annofile.open(filename);
  for(int i=0;i<numsize;i++)
    {
      if(snpgene==SNP){
	annofile>>name[i]>>chr[i]>>posstart[i];
      }else
	annofile>>name[i]>>chr[i]>>posstart[i]>>end[i];
    }
  annofile.close();
}


void initfile(const char* filename,int snptotal,int genetotal,int bsi, hsize_t* cordims,hid_t& filespace, hid_t& dset_id)
{
  hid_t file_id;

  file_id = H5Fcreate (filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  cordims = (hsize_t*)mkl_malloc(3*sizeof(hsize_t),64);

  cordims[0]=snptotal;
  cordims[1]=genetotal;
  cordims[2]=bsi;

  filespace = H5Screate_simple(3,cordims,NULL);
  
  dset_id = H5Dcreate2(file_id,"cor",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Sclose(filespace);
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
