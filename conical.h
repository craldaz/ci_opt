#ifndef CONICAL_H
#define CONICAL_H
#include "icoord.h"

class Conical {
	private:
		double* dgrad;
		double* dvec;
	
	public:
		ICoord molecule;

		//constructor to initialize molecule
		Conical(string xyzfile);
	

};

#endif
