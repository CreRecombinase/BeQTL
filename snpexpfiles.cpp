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


char* readfparam(const char*attributename,const char*attributefile)
{
  ifstream sfile;
  string ts;
  string res;
  char* cres;
  bool found=false;

  sfile.open(attributefile);
  while(getline(sfile,ts,'\n')&&(!found))
    {
      stringstream ss(ts);
      getline(ss,ts,'=');
      if(strcmp(ts.c_str(),attributename)==0)
	{
	  getline(ss,res,'=');
	  cres=(char*)malloc((res.size()+1)*sizeof(char));
	  strcpy(cres,res.c_str());
	  found=true;
	}
    }
  if(found==false)
    {
      cerr<<attributename<<" not found!!"<<endl;
    }
  return(cres);
}

    

int main(int argc, char* argv[])
{
  if(argc<2){
    cout<<"usage "<<argv[0]<<" paramfile"<<endl;
    return(-1);
  }


  double *snpmat;
  double *genemat;

  int testint;
  double testdouble;
  char* teststring;
  char* paramfile;
  
  
  string osnpcases,ogenecases,templine;
  string tsnp,tgene,tcase,ts;
  int snptotal,genetotal,casetotal;
  int snpabstotal,geneabstotal;
  int bsi,snpchunks,genechunks,cisdist;
  int snpsize,genesize;
  double t_thresh;
  char *snpfile,*genefile,*snpexph5,*annofile,*cisfilename,*transfilename,**snpnames,**genenames,*progfile;
  int i,j;

  paramfile=argv[1];
  cout<<"paramfile is "<<paramfile<<endl;

  
  snpfile=readfparam("snpfile",paramfile);

  cout<<snpfile<<endl;
  genefile=readfparam("genefile",paramfile);
  progfile=readfparam("progfile",paramfile);
  snpexph5=readfparam("h5file",paramfile);
  annofile=readfparam("annofile",paramfile);
  cisfilename=readfparam("cisfile",paramfile);
  transfilename=readfparam("transfile",paramfile);
  snpchunks=atoi(readfparam("snpchunks",paramfile));
  genechunks=atoi(readfparam("genechunks",paramfile));
  casetotal=atoi(readfparam("casetotal",paramfile));
  snptotal=atoi(readfparam("snptotal",paramfile));
  genetotal=atoi(readfparam("genetotal",paramfile));
  snpsize=atoi(readfparam("snpsize",paramfile));
  genesize=atoi(readfparam("genesize",paramfile));
  bsi=atoi(readfparam("bsi",paramfile));
  t_thresh=atof(readfparam("t_thresh",paramfile));
  cisdist=atoi(readfparam("cisdist",paramfile));
  snpabstotal=atoi(readfparam("snpabstotal",paramfile));
  geneabstotal=atoi(readfparam("geneabstotal",paramfile));
	       
  

  cout<<"parameters accepted"<<endl;


  hid_t file_id;
  hid_t snpdataspace,genedataspace,snpnamedataspace,genenamedataspace;
  hid_t snpdataset,genedataset,snpnamedataset,genenamedataset;
  hid_t fstrtype;
  hid_t root;
  herr_t status;
  hsize_t snpnamedims[1]={snptotal};
  hsize_t genenamedims[1]={genetotal};
  hsize_t snpmatsize[2]={casetotal,snptotal};
  hsize_t genematsize[2]={casetotal,genetotal};

  ifstream sfile;
  ifstream efile;


  cout<<"creating snpmat of size: "<<snptotal*casetotal<<" and expmat of size: "<<genetotal*casetotal<<endl;
  //SNPmat will NOT be transposed upon read in, but genemat will
  snpmat = (double*) malloc(snptotal*casetotal*sizeof(double));
  snpnames = (char**)malloc(snptotal*sizeof(char*));

  cout<<"reading in SNP matrix "<<snpfile<<endl;
  sfile.open(snpfile);
  getline(sfile,ts,'\n');
  i=0;
  while(getline(sfile,ts,'\n'))
    {j=0;
      stringstream ss(ts);
      getline(ss,tsnp,'\t');
      snpnames[i] = (char*)malloc((tsnp.size()+1)*sizeof(char));
      strcpy(snpnames[i],tsnp.c_str());
      while(getline(ss,tsnp,'\t'))
	{
	  snpmat[index(j,i,snptotal)]=atof(tsnp.c_str());
	  j++;
	}
      i++;
    }
  sfile.close();
  
  file_id = H5Fcreate(snpexph5,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  fstrtype= H5Tcopy(H5T_C_S1);
  H5Tset_size(fstrtype,H5T_VARIABLE);

  snpnamedataspace = H5Screate_simple(1,snpnamedims,NULL);
  snpnamedataset= H5Dcreate(file_id,"snpnames",fstrtype,snpnamedataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(snpnamedataset,fstrtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpnames);

  snpdataspace = H5Screate_simple(2,snpmatsize,NULL);
  snpdataset = H5Dcreate(file_id,"snps",H5T_NATIVE_DOUBLE,snpdataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(snpdataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpmat);


  free(snpmat);
  free(snpnames);
  cout<<"SNPmatrix Written"<<endl;

  status = H5Sclose(snpnamedataspace);
  status = H5Dclose(snpnamedataset);
  status = H5Sclose(snpdataspace);
  status = H5Dclose(snpdataset);
  
  
  cout<<"Reading Genemat"<<endl;
  genemat = (double*) malloc(casetotal*genetotal*sizeof(double));
  genenames = (char**) malloc(genetotal*sizeof(char*));
  efile.open(genefile);
  getline(efile,ts,'\n');
  i=0;
  while(getline(efile,ts,'\n'))
    {
      j=0;
      stringstream ss(ts);
      getline(ss,tgene,'\t');
      genenames[i]=(char*)malloc((tgene.size()+1)*sizeof(char));
      strcpy(genenames[i],tgene.c_str());
      while(getline(ss,tgene,'\t'))
	{

	  genemat[index(j,i,genetotal)]=atof(tgene.c_str());
	  j++;
	}
      i++;
    }
    efile.close();

    genenamedataspace = H5Screate_simple(1,genenamedims,NULL);
    genenamedataset = H5Dcreate(file_id,"genenames",fstrtype,genenamedataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(genenamedataset,fstrtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,genenames);
    genedataspace= H5Screate_simple(2,genematsize,NULL);
    genedataset = H5Dcreate(file_id,"genes",H5T_NATIVE_DOUBLE,genedataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Dwrite(genedataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,genemat);

    cout<<"Freeing genemat"<<endl;
    free(genemat);
    free(genenames);


    status = H5Sclose(genenamedataspace);
    status = H5Dclose(genenamedataset);
    status = H5Sclose(genedataspace);
    status = H5Dclose(genedataset);

    writeparameter(file_id,"annofile",annofile);
    readparam(file_id,"annofile",teststring);
    cout<<annofile<<endl;
    writeparameter(file_id,"progfile",progfile);
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
    writeparameter(file_id,"snpabstotal",snpabstotal);
    writeparameter(file_id,"geneabstotal",geneabstotal);


    writeparameter(file_id,"cisdist",cisdist);

    status = H5Fclose(file_id);
    return(0);
    
    
}
    
