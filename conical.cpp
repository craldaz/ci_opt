#include "conical.h"
using namespace std;


Conical::Conical(double* xyz,int natom,string* atoms,int* anumber,double* masses, double conv_tol, int max_steps, int run, int nprocs, int guess_wfn, string infile) 
{		
//	printf(" In the constructor\n");
	natoms=natom;
//	printf(" natoms = %i\n",natoms);
	
	anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];
	coords = new double[3*natoms+1];
	dgrad = new double[3*natoms+1];
	dvec = new double[3*natoms+1];
	
	for (int i=0;i<natoms;i++)
	{
		anames[i]=atoms[i];
		anumbers[i]=anumber[i];
		amasses[i]=masses[i];
	}
	for (int i=0;i<3*natoms;i++)
		coords[i]=xyz[i];

	CONV_TOL = conv_tol;
	STEP_OPT_ITERS=max_steps;
	runNum= run;
	ncpu=nprocs;
	molecule.alloc(natoms);
	molecule.reset(natoms,anames,anumbers,coords);
	molecule.print_xyz();
	
	create_space();

	if (guess_wfn)
		molecule.grad1.seedType = 3;
	else
  	molecule.grad1.seedType = 1; //seed from INIT1
	
  molecule.grad_init(infile,ncpu,runNum,runend,0); //level 3 is exact kNNR only, 0 is QM grad always
	V0=molecule.grad1.grads(coords, molecule.grad, molecule.Ut, 3);
	printf(" V0=%1.2f\n",V0);
  nicd0 = molecule.nicd0;
  size_ic = molecule.nbonds+molecule.nangles+molecule.ntor;
  nstates = molecule.grad1.nstates;
  wstate = molecule.grad1.wstate;
  wstate2 = molecule.grad1.wstate2;
  printf(" nstate: %i, wstate: %i, wstate2: %i\n",nstates,wstate,wstate2);

	dgradq = new double[nicd0];
	dvecq = new double[nicd0];
	dgradq_U = new double[size_ic];
	dvecq_U = new double[size_ic];

	constrain_space();
	molecule.make_Hint();
	//molecule.Hintp_to_Hint();
	printf(" finished constructing\n\n");
	
}


int Conical::create_space()
{
	int done = molecule.ic_create();
	molecule.print_ic();
	molecule.bmat_alloc();
	molecule.bmatp_create();
	molecule.bmatp_to_U();
	molecule.bmat_create();
	printf(" Done creating bmatrix stuff\n");
	return;
}

void Conical::update_space()
{
	molecule.update_ic();
	molecule.bmatp_create();
	molecule.bmatp_to_U();
	molecule.bmat_create();
#if 1
	norm_dg=molecule.dgrot_mag(dgradq,dvecq);
	molecule.project_gradq(dvecq,dvecq_U);
	molecule.project_gradq(dgradq,dgradq_U);
#else
	molecule.project_gradq(dvecq,dvecq_U);
	norm_dg=molecule.project_gradq(dgradq,dgradq_U);
	printf(" norm_dg = %1.2f",norm_dg); 
#endif

	molecule.constrain_bp(dgradq_U,dvecq_U);
	molecule.bmat_create();
	molecule.print_q();
  molecule.Hintp_to_Hint();
	return;
}

double Conical::constrain_space()
{
  for (int i=0;i<nstates-1;i++)
  {
	 molecule.grad1.dE[i] = molecule.grad1.E[i+1] - molecule.grad1.E[i];
	 printf(" dE[%i]:  %5.4f\t ",i+1,molecule.grad1.dE[i]);
  }
  printf("\n");
  dE = molecule.grad1.dE[wstate2-2]/627.5; //kcal2Hartree
#if 1
	printf(" Grads\n");
	for (int j=0;j<2;j++)
	{for (int i=0;i<3*natoms;i++)	
		printf(" %1.3f",molecule.grad1.grada[j][i]);
		printf(" \n");
	}
	
#endif
	printf(" Calculating dgrad\n");
	for (int i=0;i<3*natoms;i++)	
		dgrad[i] = molecule.grad1.grada[1][i]- molecule.grad1.grada[0][i];
	printf(" Calculating dvec[0]\n");
	molecule.grad1.dvec_calc(molecule.coords,dvec,1,0);

#if 1
	printf(" printing dgrad\n");
	for (int i=0;i<3*natoms;i++)
		printf("%1.3f\t",dgrad[i]);	
	printf("\n");
#endif
#if 0
	double test=0;
	for (int i=0;i<3*natoms;i++)
		test+=dgradi]*dgrad[i];
	test=sqrt(test);
	printf(" norm dgrad cartesians %1.2f\n",test);
#endif

	molecule.dgrad_to_q(dgrad,dgradq);
	molecule.dgrad_to_q(dvec,dvecq);

	norm_dg=0.; 
#if 0
	norm_dg=molecule.project_gradq(dgradq,dgradq_U);
	molecule.project_gradq(dvecq,dvecq_U);
#else 
	norm_dg=molecule.dgrot_mag(dgradq,dvecq);
	molecule.project_gradq(dvecq,dvecq_U);
	molecule.project_gradq(dgradq,dgradq_U);
#endif

#if 0
    printf(" printing dgradq_U \n");
	for (int j=0;j<size_ic;j++)
		printf(" %1.3f ",dgradq_U[j]);
    printf("\n");
#endif

	molecule.constrain_bp(dgradq_U,dvecq_U);
	molecule.bmat_create();
	molecule.print_q();

	return norm_dg;
}


void Conical::combined_step()
{
  for (int i=0;i<nicd0;i++)
    molecule.dq0[i] = 0.;
	update_space();
	molecule.grad_to_q();
 // molecule.dq0[nicd0-1] = -dE/norm_dg; //not sure
//  printf(" dq0[constraint]: %1.2f \n",molecule.dq0[nicd0-1]);
//	if (molecule.dq0[nicd0-1] < -0.1)
//			molecule.dq0[nicd0-1]=-0.1;
  molecule.update_ic_eigen();
  molecule.ic_to_xyz();
  molecule.print_xyz();
	exit(-1);
  energy = molecule.grad1.grads(molecule.coords, molecule.grad, molecule.Ut, 1) - V0;
	printf(" Average energy = %1.2f\n",energy);
	for (int i=0;i<nstates-1;i++)
	{
		molecule.grad1.dE[i] = molecule.grad1.E[i+1] - molecule.grad1.E[i];
		printf(" dE[%i]:  %5.4f\t ",i,molecule.grad1.dE[i]); 
	}
		printf("\n");

	return;
}

void Conical::algorithm()
{
	combined_step();
	return;
}
