#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "icoord.h"
#include "conical.h"

using namespace std;

int main(int argc, char* argv[])
{
  string inpfile;
  string run; // a number 
  string nprocs;
  switch (argc){
  case 1:
    inpfile="inpfileq";
    run="1";
    nprocs="1";
    break;
  case 2:
    inpfile="inpfileq";
    run=argv[1];
    nprocs="1";
    break;
  case 3:
    inpfile="inpfileq";
    run=argv[1];
    nprocs=argv[2];
    break;
  default:
    cout << "Invalid command line options." << endl;
    return -1;
  }

  int nnprocs = atoi(nprocs.c_str());
  printf(" Number of QC processors: %i \n",nnprocs);
  int runn = atoi(run.c_str());
	printf(" runname = %i\n",runn);
	
	
  string nstr=StringTools::int2str(runn,4,"0");
  string xyzfile = "scratch/initial"+nstr+".xyz";
  printf("  -structure filename from input: %s \n",xyzfile.c_str());


	Conical meci(xyzfile,inpfileq);
	meci.create_space();
	

	return 0;
}


