/*
	mstSmooth.c

	Copyright (C) 2012 Frederic Wilhelm <frederic.wilhelm at epfl.ch>
	Adapted from Matlab code from Benjamin Hamner <benjamin.hamner at gmail.com>
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
#include <math.h>
#include "mex.h"
#include "utils/matrix.h"
#include "utils/determineMicrostates.h"

/*************************************************************************
 *
 *		DEFINES
 *
 *************************************************************************/
#define	SQR(a)	((a)*(a))

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* \brief MEX routine for microstate smoothing
 *
 * Parameters and output descriptions concerns the contents of prhs and nlhs.
 * For more details about the prototype of the function, see the MEX-files Guide
 * of Matlab.
 *
 * \params	prhs[0]
 *		Signal (time * number of electrodes)
 *	\params	prhs[1]
 *		Microstates (number of microstates * number of electrodes)
 *	\params	prhs[2]
 *		Smoothing window b
 *	\params	prhs[3]
 *		Smoothing penalty lambda
 *	\params	prhs[4]
 *		Convergence radius
 *	\params	prhs[5]
 *		Random seed (optional)
 * 
 * 
 *	\return	plhs[0]
 *		New classification of the signal with the microstates (in 1:Nu) (time * 1)
 *	\return	plhs[1]
 *		Signal magnitude (time * 1)
 *	\return	plhs[2]
 *		Fit measure R^2
 *	\return	plhs[3]
 *		Estimated noise variance Sigmau2
 */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int t, tw, i,k;
	int Nu, Ne, Nt;
	double norm, scal, min, score;
	int Nk;
	double *signal, *microstates, *magnitude, *iS, *iM, *iMagn;
	char *signalMicrostate, *sMBuffer, *iSM, *iSMB;
	double signalNorm2, sigmad2, sigmau2, e, R2, sigma02=0, deviation;
	double epsilon, lambda; /* Radius of convergence and smoothing penalty */
	int b; /* Smoothing window */
	bool stop = false; /* Stop the loop */

	if (nrhs < 5)
		mexErrMsgTxt("Must have at least four input argument");		

	/* Retrieve the seed */
	if (nrhs >= 6)
		srand((int) *mxGetPr(prhs[5]));
	else
		srand(time(0));

	/* Retrieve the signal data AND TRANSPOSE IT (for later convenience) */
	Nt = mxGetM(prhs[0]);
	Ne = mxGetN(prhs[0]);
	signal = malloc(Nt * Ne * sizeof(double));
	memcpy(signal, mxGetPr(prhs[0]), Ne*Nt*sizeof(double));
	transpose(signal, Nt, Ne);

	/* Retrieve the microstates */
	Nu = mxGetN(prhs[1]);
	microstates = malloc(Nu * Ne * sizeof(double));
	memcpy(microstates, mxGetPr(prhs[1]), Nu*Ne*sizeof(double));

	/* Retrieve b, lambda and the radius of convergence */
	b = (int) *mxGetPr(prhs[2]);
	lambda = *mxGetPr(prhs[3]);
	epsilon = *mxGetPr(prhs[4]);

	/* Determine the initial microstate estimate */
	signalMicrostate = malloc(Nt*sizeof(char));
	sMBuffer = malloc(Nt*sizeof(char));
	magnitude = malloc(Nt*sizeof(double));
	determineMicrostates(Nu, Ne, Nt, signal, microstates, signalMicrostate, magnitude);

	/* Determine the sum of the norms of the datas and sigmad */
	iS = signal;
	signalNorm2 = 0;
	for (k = 0 ; k < Ne*Nt ; k++)
	{
		signalNorm2 += SQR(*iS);
		iS++;
	}
	sigmad2 = signalNorm2 / (Nt * (Ne - 1));

	/* Compute sigmau */
	sigmau2 = signalNorm2;
	iMagn = magnitude;
	for (t = 0 ; t < Nt ; t++)
	{
		sigmau2 -= SQR(*iMagn);
		iMagn++;
	}
	sigmau2 = sigmau2 / (Nt * (Ne - 1));
	e = sigmau2;

	/* Iterate until converged */
	while (!stop)
	{
		/* Loop over time (eliminate border where one can't define
		the window) */
		iS = signal + b*Ne;
		iSMB = sMBuffer + b;
		for (t = b ; t < Nt-b ; t++)
		{
			/* Norm of state t */
			norm = 0;
			for (i = 0 ; i < Ne ; i++)
			{
				norm += SQR(*iS);
				iS++;
			}

			/* Examine each state */
			min = FLT_MAX;
			iM = microstates;
			for (k = 0 ; k < Nu ; k++)
			{
				/* Count the number of states k in the window */
				Nk = 0;
				for (tw = t-b ; tw <= t+b ; tw++)
					if (signalMicrostate[tw] == k)
						Nk++;
				
				/* Compute score of state k for time t */
				score = norm;
				iS -= Ne;
				scal = 0;
				for (i = 0 ; i < Ne ; i++)
				{
					scal += (*iM++)*(*iS++);
				}
				score = score - SQR(scal);
				score = (score / (2*e*(Ne-1))) - lambda*Nk;
				
				/* Calculate min */
				if (score < min)
				{
					min = score;
					*iSMB = k;
				}
			}

			/* Next time */
			iSMB++;
		}

		/* Update signalMicrostate */
		memcpy(signalMicrostate+b, sMBuffer+b, (Nt-2*b)*sizeof(char));

		/* Compute sigmau */
		sigmau2 = signalNorm2;
		iS = signal;
		iSM = signalMicrostate;
		for (t = 0 ; t < Nt ; t++)
		{
			iM = microstates + (*iSM)*Ne;
			scal = 0;
			for (i = 0 ; i < Ne ; i++)
				scal += (*iM++)*(*iS++);
			sigmau2 -= SQR(scal);
			iSM++;
		}
		sigmau2 = sigmau2 / (Nt * (Ne - 1));

		/* Compute deviation and R2 */
		deviation = fabs(sigma02 - sigmau2) / sigmau2;
		R2 = 1 - sigmau2/sigmad2;

		/* Examine convergence */
		if (deviation < epsilon)
			stop = true;
		else
			sigma02 = sigmau2;
	}

	/* Create the output for the classification of the signal */
	if (nlhs >= 1)
	{
		iSM = signalMicrostate;
		for (t = 0 ; t < Nt ; t++)
			(*iSM++) += 1;
		plhs[0] = mxCreateNumericMatrix(Nt, 1, mxUINT8_CLASS, mxREAL);
		memcpy(mxGetPr(plhs[0]), signalMicrostate, Nt*sizeof(char));
	}

	/* Create the output for the magnitude (which has to be computed) */
	if (nlhs >= 2)
	{
		iMagn = magnitude;
		iS = signal;
		iSM = signalMicrostate;
		for (t = 0 ; t < Nt ; t++)
		{
			scal = 0;
			iM = microstates + (*iSM++ - 1)*Ne; /* <!> previous conversion to
															 Matlab's convention */
			for (i = 0 ; i < Ne ; i++)
				scal += (*iM++) * (*iS++);
			*iMagn++ = scal;
		}
		plhs[1] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
		memcpy(mxGetPr(plhs[1]), magnitude, Nt*sizeof(double));
	}

	/* Create the output for the measure R2 */
	if (nlhs >= 3)
		plhs[2] = mxCreateDoubleScalar(R2);

	/* Create the output for the measure R2 */
	if (nlhs >= 4)
		plhs[3] = mxCreateDoubleScalar(sigmau2);

	/* Relesase memory */
	free(signal);
	free(microstates);
	free(magnitude);
	free(signalMicrostate);
	free(sMBuffer);
}
