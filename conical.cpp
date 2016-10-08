#include "conical.h"
using namespace std;


//Conical::Conical( string xyzfile, string infilename)
//{
//  printf("Initializing Tolerances and Parameters... \n");
//  printf("  -Opening %s \n",infilename.c_str());
//
//  ifstream infile;
//  infile.open(infilename.c_str());
//  if (!infile){
//    printf("\n Error opening %s \n",infilename.c_str());
//    exit(-1);
//  }
//
//  printf("  -reading file... \n");
//  // pass infile to stringtools and get the line containing tag
//  string tag="String Info";
//  bool found=StringTools::findstr(infile, tag);
//  if (!found) { cout << "Could not find tag for Default Info" << endl; exit(-1);}
//  string line, templine, tagname;
//
//  // parse the input section
//  bool stillreading=true;
//  while (stillreading)
//  {
//    stillreading=false;
//    // set to false and set back to true if we read something
//    // get filename
//    getline(infile, line);
//    vector<string> tok_line = StringTools::tokenize(line, " ,\t");
//    templine=StringTools::newCleanString(tok_line[0]);
//    tagname=StringTools::trimRight(templine);
//    // these variables are denoted by strings with same name
//    if (tagname=="STEP_OPT_ITERS") {
//      STEP_OPT_ITERS = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -STEP_OPT_ITERS: " << STEP_OPT_ITERS << endl;
//    }
//    if (tagname=="CONV_TOL") {
//      CONV_TOL=atof(tok_line[1].c_str());
//      stillreading=true;
//      cout <<"  -CONV_TOL = " << CONV_TOL << endl;
//    }
//    if (tagname=="ADD_NODE_TOL"){
//      ADD_NODE_TOL=atof(tok_line[1].c_str());
//      stillreading=true;
//      cout <<"  -ADD_NODE_TOL = " << ADD_NODE_TOL << endl;
//    }
//	  if (tagname=="RESTART_WFN"){
//			restart_wfn = atoi(tok_line[1].c_str());
//			stillreading=true;
//			cout <<"  -restart_wfn = " << restart_wfn << endl; 
//		}
//    if (tagname=="nnodes" || tagname=="NNODES"){
//      nnmax = atoi(tok_line[1].c_str());
//      stillreading = true;
//      cout <<"  -NNODES = " << nnmax << endl;
//      if (nnmax < 3)
//      {
//        printf("\n\n ERROR: NNODES cannot be less than 3 \n");
//        exit(1);
//      }
//    }
//
//  } //while stillreading
//  infile.close();
//  nnmax0 = nnmax;
//
//  printf(" Done reading inpfileq \n\n");
//	molecule.init(xyzfile);
//	molecule.grad_init(infileq,ncpu,run,0,0);
//}


int Conical::create_space()
{
	int done = molecule.ic_create();
	molecule.print_ic();
	molecule.bmat_alloc();
	molecule.bmatp_create();
	molecule.bmatp_to_U();
	molecule.bmat_create();
	printf(" Done creating bmatrix stuff\n");

}
