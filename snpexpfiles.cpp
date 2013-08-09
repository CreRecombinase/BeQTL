#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <sstream>
#include <set>
#include <map>
#include <string>
using namespace std;

int index (int i, int j, int M)
{
  return((M*i)+j);
}


void writeparameter(hid_t file_id,const char* attributename, const char* attributevalue){
  hid_t root,dataspace,att,fstrtype;
  hsize_t dim[1]={1};
  herr_t ret;

  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  dataspace = H5Screate_simple(1,dim,NULL);
  fstrtype= H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);
  att=H5Acreate(root,attributename,fstrtype,dataspace,H5P_DEFAULT,H5P_DEFAULT);
  ret=H5Awrite(att,fstrtype,&attributevalue);
  ret=H5Sclose(dataspace);
  ret=H5Aclose(att);
  


  
  
}

void writeparameter(hid_t file_id,const char* attributename, const int attributevalue){
  hid_t root,dataspace,att;
  hsize_t dim[1]={1};
  herr_t ret;
  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  dataspace=H5Screate_simple(1,dim,NULL);
  att=H5Acreate(root,attributename,H5T_NATIVE_INT,dataspace,H5P_DEFAULT,H5P_DEFAULT);
  ret=H5Awrite(att,H5T_NATIVE_INT,&attributevalue);
  ret=H5Sclose(dataspace);
  ret=H5Aclose(att);
  
}

void writeparameter(hid_t file_id,const char* attributename, const double attributevalue){
  hid_t root,dataspace,att;
  hsize_t dim[1]={1};
  herr_t ret;
  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  dataspace=H5Screate_simple(1,dim,NULL);
  att=H5Acreate(root,attributename,H5T_NATIVE_DOUBLE,dataspace,H5P_DEFAULT,H5P_DEFAULT);
  ret=H5Awrite(att,H5T_NATIVE_DOUBLE,&attributevalue);
  ret=H5Sclose(dataspace);
  ret=H5Aclose(att);
  
}


void readparam(hid_t file_id, const char* attributename,char*&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;

  root=H5Gopen(file_id,"/",H5P_DEFAULT);
  fstrtype= H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,fstrtype,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}
	  
void readparam(hid_t file_id,const char* attributename,int&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,H5T_NATIVE_INT,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}

void readparam(hid_t file_id,const char* attributename,double&attributevalue)
{
  hid_t root,att;
  hid_t fstrtype;
  
  root = H5Gopen(file_id,"/",H5P_DEFAULT);
  att=H5Aopen_name(root,attributename);
  H5Aread(att,H5T_NATIVE_DOUBLE,&attributevalue);
  H5Gclose(root);
  H5Aclose(att);
}

int main(int argc, char* argv[])
{
  if(argc<4){
    cout<<"usage "<<argv[0]<<" snpfile genefile h5file annofile cisfilename transfilename snpchunks genechunks casetotal snptotal genetotal snpsize genesize bsi t_thresh cisdist"<<endl;
    return(-1);
  }

  char **snpnames,**genenames;

  double *snpmat;
  double *genemat;

  int testint;
  double testdouble;
  char* teststring;
  
  string osnpcases,ogenecases,templine;
  string tsnp,tgene,tcase,ts;
  int snptotal,genetotal,casetotal;
  int bsi,snpchunks,genechunks,cisdist;
  int snpsize,genesize;
  double t_thresh;
  char *snpfile,*genefile,*snpexph5,*annofile,*cisfilename,*transfilename,*readfilelock,*writefilelock;
  int i,j;

  
  snpfile=argv[1];
  genefile=argv[2];
  snpexph5=argv[3];
  annofile=argv[4];
  cisfilename=argv[5];
  transfilename=argv[6];
  snpchunks=atoi(argv[7]);
  genechunks=atoi(argv[8]);
  casetotal=atoi(argv[9]);
  snptotal=atoi(argv[10]);
  genetotal=atoi(argv[11]);
  snpsize=atoi(argv[12]);
  genesize=atoi(argv[13]);
  bsi=atoi(argv[14]);
  t_thresh = atof(argv[15]);
  cisdist=atoi(argv[16]);
  

  cout<<"parameters accepted"<<endl;


  hid_t file_id;
  hid_t snpdataspace,genedataspace;
  hid_t snpdataset,genedataset;
  hid_t fstrtype;
  hid_t root;
  herr_t status;
  hsize_t snpmatsize[2]={casetotal,snptotal};
  hsize_t genematsize[2]={casetotal,genetotal};

  ifstream sfile;
  ifstream efile;
  istringstream ssnpc,sgenec,cases;


  cout<<"creating snpmat of size: "<<snptotal*casetotal<<" and expmat of size: "<<genetotal*casetotal<<endl;
  //SNPmat will NOT be transposed upon read in, but genemat will
  snpmat = (double*) malloc(snptotal*casetotal*sizeof(double));

  cout<<"reading in SNP matrix "<<snpfile<<endl;
  sfile.open(snpfile);
  getline(sfile,ts,'\n');
  i=0;
  while(getline(sfile,ts,'\n'))
    {j=0;
      stringstream ss(ts);
      getline(ss,tsnp,'\t');
      while(getline(ss,tsnp,'\t'))
	{
	  snpmat[index(j,i,snptotal)]=atof(tsnp.c_str());
	  j++;
	}
      i++;
    }
  sfile.close();
  
  file_id = H5Fcreate(snpexph5,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  snpdataspace = H5Screate_simple(2,snpmatsize,NULL);
  snpdataset = H5Dcreate(file_id,"snps",H5T_NATIVE_DOUBLE,snpdataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(snpdataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpmat);


  free(snpmat);
  cout<<"SNPmatrix Written"<<endl;

  status = H5Sclose(snpdataspace);
  status = H5Dclose(snpdataset);
  
  
  cout<<"Reading Genemat"<<endl;
  genemat = (double*) malloc(casetotal*genetotal*sizeof(double));
  efile.open(genefile);
  getline(efile,ts,'\n');
  i=0;
  while(getline(efile,ts,'\n'))
    {
      j=0;
      stringstream ss(ts);
      getline(ss,tgene,'\t');
      while(getline(ss,tgene,'\t'))
	{

	  genemat[index(j,i,genetotal)]=atof(tgene.c_str());
	  j++;
	}
      i++;
    }
    efile.close();

    genedataspace= H5Screate_simple(2,genematsize,NULL);
    genedataset = H5Dcreate(file_id,"genes",H5T_NATIVE_DOUBLE,genedataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(genedataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,genemat);

    cout<<"Freeing genemat"<<endl;
    free(genemat);
    
    status = H5Sclose(genedataspace);
    status = H5Dclose(genedataset);

    writeparameter(file_id,"annofile",annofile);
    readparam(file_id,"annofile",teststring);
    cout<<annofile<<endl;
    writeparameter(file_id,"cisfilename",cisfilename);
    writeparameter(file_id,"transfilename",transfilename);
    writeparameter(file_id,"snpchunks",snpchunks);
    writeparameter(file_id,"genechunks",genechunks);
    writeparameter(file_id,"casetotal",casetotal);
    writeparameter(file_id,"snptotal",snptotal);
    writeparameter(file_id,"genetotal",genetotal);
    writeparameter(file_id,"bsi",bsi);
    writeparameter(file_id,"t_thresh",t_thresh);
    writeparameter(file_id,"snpsize",snpsize);
    writeparameter(file_id,"genesize",genesize);

    writeparameter(file_id,"cisdist",cisdist);

    status = H5Fclose(file_id);
    return(0);
    
    
}
    
