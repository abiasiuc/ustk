/*
	mstGaussian.c

	Copyright (C) 2010-2015 
	Andrea Biasiucci <andrea.biasiucci at epfl.ch>	
	Frederic Wilhelm <frederic.wilhelm at epfl.ch>
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
#include "utils/gaussian.h"

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* \brief MEX routine accounting for a gaussian (microstate model)
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
 *		Microstate (center of gaussian) (Ns)
 *	\params	prhs[2]
 *		Variance (double)
 * 
 *	\return	plhs[0]
 *		Responsibilities (Nt * Nu)
  */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int Nu, Ns, Nt, t;
	double *signal, *iS,
			 *microstate,
			 sigma2,
			 *probabilities, *iP;
	
	/* Retrieve the signal data AND TRANSPOSE IT (for later convenience) */
	Nt = mxGetM(prhs[0]);
	Ns = mxGetN(prhs[0]);
	signal = malloc(Ns*Nt*sizeof(double));
	memcpy(signal, mxGetPr(prhs[0]), Ns*Nt*sizeof(double));
	transpose(signal, Nt, Ns);

	/* Retrieve the microstate */
	microstate = mxGetPr(prhs[1]);

	/* Retrieve the variance */
	sigma2 = *mxGetPr(prhs[2]);

	/* Compute the probabilities */
	if (nlhs >= 1)
	{
		/* Create the output (probabilities) */
		plhs[0] = mxCreateDoubleMatrix(Nt, 1, mxREAL);

		/* Compute the probabilities */
		probabilities = mxGetPr(plhs[0]);
		iP = probabilities;
		iS = signal;
		for (t=0 ; t<Nt ; t++)
		{
			*iP++ = gaussian(iS, microstate, sigma2, Ns);
	
			/* Next time sample */
			iS += Ns;
		}
	}

	/* Release memory */
	free(signal);
}
