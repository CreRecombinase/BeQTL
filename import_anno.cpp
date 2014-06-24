//Code to import annotation data
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;

int index (int i, int j, int M)
{
  return((M*i)+j);
}

int main(int argc, char* argv[])
{

  char* snpannofile,geneannofile,outsnpfile;
  int snpsize,genesize;
  snpannofile = argv[1];
  geneannofile = argv[2];
  outsnpfile = argv[3];
  genesize = atoi(argv[4]);
  snpsize = atoi(argv[5]);


  char snpannofile[]="/home/nwk2/mkl_test/snpanno.txt";
  char geneannofile[]="/home/nwk2/mkl_test/geneanno.txt";

  char outsnpfile[]="/home/nwk2/mkl_test/snpgeneanno.h5";

  
  hsize_t snpnamedims[1]={snpsize};
  hsize_t genenamedims[1]={genesize};
  hid_t fstrtype;
  hid_t mstrtype;
  hid_t file_id;
  hid_t snpnamespace,snparrayspace;
  hid_t genenamespace,genearrayspace;
  hid_t snpnamedset,snparraydset;
  hid_t genenamedset,genearraydset;
  hsize_t snpdims[2],genedims[2];
  herr_t status;

  char **snpname;
  char **genename;
  int* snparray;
  int* genearray;
  ifstream annofile;
  string tstring;

  int snp_colsize =2;
  int gene_colsize=3;

  snparray = (int *)malloc(2*snpsize*sizeof(int));
  snpname = (char**)malloc(snpsize*sizeof(char*));
  
  


  annofile.open(snpannofile);
  for(int i=0;i<snpsize;i++)
    {
      annofile>>tstring>>snparray[index(i,0,snp_colsize)]>>snparray[index(i,1,snp_colsize)];
      snpname[i] = (char*)malloc(tstring.length()+1*sizeof(char));
      strcpy(snpname[i],tstring.c_str());
      //cout<<snpname[i]<<"\t"<<snparray[index(i,0,snpsize)]<<"\t"<<snparray[index(i,1,snpsize)]<<endl;
    }
  annofile.close();



  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  mstrtype = H5Tcopy(H5T_C_S1);
  status= H5Tset_size(mstrtype,H5T_VARIABLE);

  
  
  
  file_id = H5Fcreate(outsnpfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  snpnamespace = H5Screate_simple(1,snpnamedims,NULL);
  
  snpnamedset = H5Dcreate(file_id,"snpname",fstrtype,snpnamespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  status=H5Dwrite(snpnamedset,mstrtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpname);
  //close spaces create for snpname array
  status = H5Dclose(snpnamedset);
  status= H5Sclose(snpnamespace);
  free(snpname);
  
  snpdims[0]=snpsize;
  snpdims[1]=2;
  
  snparrayspace = H5Screate_simple(2,snpdims,NULL);

  snparraydset = H5Dcreate(file_id,"snparray",H5T_NATIVE_INT,snparrayspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  status=H5Dwrite(snparraydset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,snparray);
  //close snparray spaces and close snpannofile
  status = H5Dclose(snparraydset);
  status=H5Sclose(snparrayspace);
  free(snparray);

  
  //start reading in and writing gene annofile
  genedims[0]=genesize;
  genedims[1]=3;

  genearray=(int *)malloc(3*genesize*sizeof(int));
  genename=(char**)malloc(genesize*sizeof(char*));

  annofile.open(geneannofile);
  for(int i=0;i<genesize;i++)
    {
      annofile>>tstring>>genearray[index(i,0,gene_colsize)]>>genearray[index(i,1,gene_colsize)]>>genearray[index(i,2,gene_colsize)];
      genename[i]=(char*)malloc((tstring.length()+1)*sizeof(char));
      strcpy(genename[i],tstring.c_str());
      //cout<<genename[i]<<"\t"<<genearray[index(i,0,gene_colsize)]<<"\t"<<genearray[index(i,1,gene_colsize)]<<"\t"<<genearray[index(i,2,gene_colsize)]<<endl;
    }
  annofile.close();
  
 
  
  genenamespace = H5Screate_simple(1,genenamedims,NULL);
  genenamedset = H5Dcreate(file_id,"genename",fstrtype,genenamespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(genenamedset,mstrtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,genename);
  status=H5Dclose(snpnamedset);
  status=H5Sclose(snpnamespace);
  free(genename);

  genearrayspace = H5Screate_simple(2,genedims,NULL);
  genearraydset = H5Dcreate(file_id,"genearray",H5T_NATIVE_INT,genearrayspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(genearraydset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,genearray);

  status = H5Tclose(fstrtype);
  status=H5Tclose(mstrtype);
  status= H5Dclose(genearraydset);
  status=H5Sclose(genearrayspace);
  status=H5Fclose(file_id);
  free(genearray);

  return(0);
}
