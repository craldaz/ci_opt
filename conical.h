#ifndef CONICAL_H
#define CONICAL_H
#include "icoord.h"

class Conical {
	private:
		double* dgrad;
		double* dvec;
		double* dgradq;
		double* dvecq;
		double* dgradq_U;
		double* dvecq_U;
		double dE; // energy gap in kcal/mol
		double* coords;
		int natoms;
  	int ncpu; // number of cpu for qchem,molpro, etc
  	int runNum; // unique id for this gstring instance
		int runend; 
  	int STEP_OPT_ITERS;
  	double CONV_TOL;
  	string* anames;		//array of atomic symbols (for creating input QC file)
  	int* anumbers;		//array of atomic indices (for looking up period table stuff)
  	double* amasses;		//array of atomic masses (used for mass-weighting coordinates)
		double V0;			// Average Energy at beginning
		
		int size_ic; // number of primitive coordinates
		int nicd0; // number of delocalized internal coordinates
		int nstates;
		int wstate;
		int wstate2;
		double energy;
		double norm_dg;
		
	
	public:
		ICoord molecule;
		//constructor to initialize molecule
		 Conical(double* xyz,int natoms,string* atoms, int* anumber, double* masses, double conv_tol, int max_steps, int run, int ncpu,int guess_wfn,string infile);
		 //Conical();
		int create_space();
		double constrain_space();
		void update_space();
		void combined_step();
		void algorithm();
	

};

#endif
