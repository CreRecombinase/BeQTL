//Code to import annotation data

#include <ifstream>
#include <sstream>
#include <string>
#include <string.h>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;

int index (int i, int j, int M)
{
  return((M*i)+j);
}

int main()
{
  char snpannofile[]="/home/nwk2/mkl_test/snpanno.txt";
  char geneannofile[]="/home/nwk2/mkl_test/geneanno.txt";

  char outsnpfile[]="/home/nwk2/mkl_test/snpanno.h5";
  char outgenefile[]="/home/nwk2/mkl_test/geneanno.h5";
  
  int snpsize[1]={906598};
  int genesize[1]={20501};
  hid_t fstrtype;
  hid_t mstrmtype;
  hid_t snpfile_id,genefile_id;
  hid_t snpnamespace,snparrayspace;
  hid_t snpnamedset,snparraydset;
  hsize_t snpdims[2],genedims[2];
  herr_t status;

  char **snpname;
  int* snpchr,*snppos;
  int* genechr,*genestart,*geneend;
  ifstream annofile;
  strting tstring;

  snparray = (int *)malloc(2*snpsize[0]*sizeof(int));
  snpname = (char**)malloc(snpsize[0]*sizeof(char*));
  
  


  annofile.open(snpannofile);
  for(int i=0;i<snpsize[0];i++)
    {
      annofile>>tstring>>snparray[index(i,0,snpsize[0])]>>snparray[index(i,1,snpsize)];
      snpname[i] = (char*)malloc(tstring.length()+1*sizeof(char));
      strcpy(snpname[i],tstring.c_str());
    }
  annofile.close();



  fstrtype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(fstrtype,H5T_VARIABLE);
  mstrtype = H5Tcopy(H5T_C_S1);
  status= H5Tset_size(mstrtype,H5T_VARIABLE);

  
  
  
  snpfile_id = H5Fcreate(outsnpfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  snpnamespace = H5Screate_simple(1,snpsize,NULL);
  
  snpnamedset = H5Dcreate(snpfile_id,"snpname",fstrtype,snpnamespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  


  

  




  genefile_id = H5Fcreate(geneannofile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  
