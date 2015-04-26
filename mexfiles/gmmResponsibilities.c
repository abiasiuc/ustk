/*
	gmmResponsibilities.c

	Copyright (C) 2010-2015
	Andrea Biasiucci <andrea.biasiucci at epfl.ch> 
	Frederic Wilhelm <frederic.wilhelm at epfl.ch>
	Benjamin Hamner <benjamin.hamner at gmail.com>
	Chair in Non-Invasive Brain-Machine Interface (http://cnbi.epfl.ch/)
	EPFL Ecole Polytechnique Federale de Lausanne (http://www.epfl.ch)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	*/

#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "utils/matrix.h"
#include "utils/determineMicrostates.h"

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* \brief MEX routine for computing the responsibilities of signals
 *
 * Parameters and output descriptions concerns the contents of prhs and nlhs.
 * For more details about the prototype of the function, see the MEX-files Guide
 * of Matlab.
 *
 * Notations:	Nt: Number of time samples
 *					Ns: Number of electrodes
 *					Nu: Number of microstates
 *
 * \params	prhs[0]
 *		Signal (Nt * Ns)
 *	\params	prhs[1]
 *		Microstates (center of gaussians) (Ns * Nu)
 *	\params	prhs[2]
 *		Variances (1 * Nu)
 *  \params prhs[3]
 *    Mixing coefficients (1 * Nu)
 * 
 *	\return	plhs[0]
 *		Responsibilities (Nt * Nu)
  */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int Nu, Ns, Nt;
	double *signal,
			 *microstates,
			 *sigma2,
			 *mix,
			 *gamma;

	if (nrhs < 4)
		mexErrMsgTxt("Must have at least four input argument");

	/* Retrieve the dimensions */
	Nt = mxGetM(prhs[0]);
	Ns = mxGetN(prhs[0]);
	Nu = mxGetN(prhs[1]);

	/* Retrieve the datas */
	signal = malloc(Ns*Nt*sizeof(double));
	memcpy(signal, mxGetPr(prhs[0]), Ns*Nt*sizeof(double));
	transpose(signal, Nt, Ns); /* For latter convenience */
	microstates = mxGetPr(prhs[1]);
	sigma2 = mxGetPr(prhs[2]);
	mix = mxGetPr(prhs[3]);

	/* Compute the responsibilities */
	if (nlhs >= 1)
	{
		/* Compute gamma */
		gamma = malloc(Nu*Nt*sizeof(double));
		computeResponsibilities(Nu, Ns, Nt, signal, microstates, sigma2,
				mix, gamma);
		transpose(gamma, Nu, Nt);

		/* Create the output */
		plhs[0] = mxCreateDoubleMatrix(Nt, Nu, mxREAL);
		memcpy(mxGetPr(plhs[0]), gamma, Nt*Nu*sizeof(double));

		/* Release gamma */
		free(gamma);
	}

	/* Release memory */
	free(signal);
}

