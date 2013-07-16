#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <sstream>

using namespace std;

int index (int i, int j, int M)
{
  return((M*i)+j);
}

int main(int argc, char* argv[])
{
  double *snpmat;
  double *expmat;
  string line,ts;
  int rows,genes,snps;
  char snpname[]="/scratch/nwk2/hdf5/snpmat_BRCAN.txt";
  char expname[]="/scratch/nwk2/hdf5/expmat_BRCAN.txt";

  char snpexph5[]="/scratch/nwk2/hdf5/snpexpmat_BRCAN.h5";
  rows = 177;
  snps = 906598;
  genes = 20501;




  snpmat = (double *) malloc(rows*snps*sizeof(double));
  expmat = (double *) malloc(rows*genes*sizeof(double));

  ifstream sfile;
  ifstream efile;
  
  sfile.open(snpname);
  efile.open(expname);

  int i=0;
  int j=0;
  cout<<"Reading in SNP data"<<endl;
  getline(sfile,line);
  while(getline(sfile,line))
    {
      j=0;
      stringstream ss(line);
      while(getline(ss,ts,'\t'))
	{
	  snpmat[index(i,j,snps)]=atof(ts.c_str());
	  j++;
	}
      i++;
    }
  sfile.close();
  cout<<"Reading in expfile"<<endl;
  i=0;
  getline(efile,line);
  while(getline(efile,line))
    {
      j=0;
      stringstream ss(line);
      while(getline(ss,ts,'\t'))
	{
	  expmat[index(i,j,genes)]=atof(ts.c_str());
	  j++;
	}
      i++;
    }
  efile.close();

  hsize_t snpmatsize[2]={rows,snps};
  hsize_t expmatsize[2]={rows,genes};
  hid_t file_id;
  hid_t snpdataspace;
  hid_t expdataspace;
  hid_t snpdataset,expdataset;
  herr_t status;

  file_id = H5Fcreate (snpexph5,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);



  snpdataspace = H5Screate_simple(2,snpmatsize,NULL);
  expdataspace = H5Screate_simple(2,expmatsize,NULL);


  snpdataset = H5Dcreate2(file_id,"snps",H5T_NATIVE_DOUBLE,snpdataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  expdataset = H5Dcreate2(file_id,"genes",H5T_NATIVE_DOUBLE,expdataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  

  status = H5Dwrite( snpdataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpmat);
  status = H5Dwrite( expdataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,expmat);



  H5Sclose(snpdataspace);
  H5Sclose(expdataspace);

  H5Dclose(snpdataset);
  H5Dclose(expdataset);
  H5Fclose(file_id);

  free(snpmat);
  free(expmat);


  return(0);
}
