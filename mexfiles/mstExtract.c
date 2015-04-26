/*
	mstExtractMicrostate.c

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

/* \brief MEX routine for microstate extraction
 *
 * Parameters and output descriptions concerns the contents of prhs and nlhs.
 * For more details about the prototype of the function, see the MEX-files Guide
 * of Matlab.
 *
 * \params	prhs[0]
 *		Signal (time * number of electrodes)
 *	\params	prhs[1]
 *		Number of states Nu (i.e. clusters)
 *	\params	prhs[2]
 *		Radius of convergence
 *	\params	prhs[3]
 *		(optional) Random seed
 *	\params	prhs[4]
 *		(optional) Initial microstates
 * 
 *	\return	plhs[0]
 *		List of states (number of electrodes * number of states)
 *	\return	plhs[1]
 *		Classification of the signal with the microstates (in 1:Nu) (time * 1)
 *	\return	plhs[2]
 *		Signal magnitude (time * 1)
 *	\return	plhs[3]
 *		Fit measure R^2
 *	\return	plhs[4]
 *		Estimated noise variance Sigmau2
 */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int t, i,k;
	int Nu, Ne, Nt;
	double norm, scal;
	double *signal, *microstates, *iS, *iS2, *iM;
	char *signalMicrostate, *iSM;
	double *magnitude, *iMagn;
	double *matS, *iMatS;
	int *numSignalInState;
	double signalNorm2, sigmad2, sigmau2, R2, sigma02=0, deviation;
	double epsilon; /* Radius of convergence */
	bool stop=false; /* stop the loop */

	if (nrhs < 3)
		mexErrMsgTxt("Must have at least three input argument");		

	/* Retrieve the seed */
	if (nrhs >= 4)
		srand((int) *mxGetPr(prhs[3]));
	else
		srand(time(0));

	/* Retrieve the signal data AND TRANSPOSE IT (for later convenience) */
	Nt = mxGetM(prhs[0]);
	Ne = mxGetN(prhs[0]);
	signal = malloc(Nt * Ne * sizeof(double));
	memcpy(signal, mxGetPr(prhs[0]), Ne*Nt*sizeof(double));
	transpose(signal, Nt, Ne);

	/* Retrieve the number of microstates */
	Nu = (int) *mxGetPr(prhs[1]);

	/* Retrieve the radius of convergence */
	epsilon = *mxGetPr(prhs[2]);
	
	/* Initialize the microstates */
	microstates = malloc(Nu * Ne * sizeof(double));
	if (nrhs >= 5)
		memcpy(microstates, mxGetPr(prhs[4]), Nu*Ne*sizeof(double));
	else
	{
		/* Initialize iterator */
		iM = microstates;
		for (i = 0 ; i < Nu ; i++)
		{
			/* Generate one random microstate */
			scal = 0; /* Average (for transformation to ARD) */
			norm = 0; /* Norm (for normalization) */
			for (k = 0 ; k < Ne ; k++)
			{
				*iM = rand() % 1000;
				scal += *iM;
				iM++;
			}
			/* Transform to average reference data */
			iM -= Ne;
			scal /= Ne;
			for (k = 0 ; k < Ne ; k++)
			{
				*iM -= scal;
				norm += SQR(*iM);
				iM++;
			}
			/* Normalize the microstate */
			iM -= Ne;
			norm = sqrt(norm);
			for (k = 0 ; k < Ne ; k++)
			{
				*iM++ /= norm;
			}
		}
	}

	/* Determine the microstates of the data */
	signalMicrostate = malloc(Nt*sizeof(char));
	magnitude = malloc(Nt*sizeof(double));
	determineMicrostates(Nu, Ne, Nt, signal, microstates, signalMicrostate, magnitude);

	/* Determine the sum of the norms of the datas and sigmad2 */
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

	/* Compute R^2 */
	R2 = 1 - sigmau2/sigmad2;

	/* Iterate until converged */
	matS = malloc(Ne*Ne*Nu*sizeof(double)); /* Matrix S (for each state) */
	numSignalInState = malloc(Nu*sizeof(int));
	while (!stop)
	{
		/* Reinitialize matrix S and numSignalInState; */
		memset(matS, 0, Ne*Ne*Nu*sizeof(double));
		memset(numSignalInState, 0, Nu*sizeof(int));

		/* Compute matrices S and numSignalInState */
		iS = signal;
		iSM = signalMicrostate;
		for (t = 0 ; t < Nt ; t++)
		{
			/* Add signal at t to corresponding S */
			iMatS = matS + (*iSM)*Ne*Ne;
			iS2 = iS;
			for (i = 0 ; i < Ne ; i++)
			{
				for (k = 0 ; k < Ne ; k++)
					(*iMatS++) += (*iS) * (*iS2++);
				iS++;
				iS2 -= Ne;
			}

			/* Increment numSignalInState */
			numSignalInState[*iSM]++;

			/* Jump to next time */
			iSM++;
		}

		/* Moves or create each microstate */
		for (k = 0 ; k < Nu ; k++)
		{
			if (numSignalInState[k] == 0)
			{
				/* Reset the microstate to a random vector */
				norm = 0; /* Norm (for normalization) */
				scal = 0; /* Average (for ARD transformation) */
				iM = microstates + k*Ne;
				for (i = 0 ; i < Ne ; i++)
				{
					*iM = rand() % 1000;
					scal += *iM;
					iM++;
				}
				/* Transform to avg reference data */
				iM -= Ne;
				scal /= Ne;
				for (i = 0 ; i < Ne ; i++)
				{
					*iM -= scal;
					norm += SQR(*iM);
					iM++;
				}
				/* Normalize microstate */
				iM -= Ne;
				norm = sqrt(norm);
				for (i = 0 ; i < Ne ; i++)
				{
					*iM++ /= norm;
				}
			}
			else
			{
				/* Computes eigenvector of largest eigenvalue of matrix S */
				iS = matS + k*Ne*Ne;
				largestEigenvector(iS, Ne, microstates + k*Ne);
			}
		}
		
		/* Update signal microstates */
		determineMicrostates(Nu, Ne, Nt, signal, microstates, signalMicrostate, magnitude);
		
		/* Compute sigmau */
		sigmau2 = signalNorm2;
		iMagn = magnitude;
		for (t = 0 ; t < Nt ; t++)
		{
			sigmau2 -= SQR(*iMagn);
			iMagn++;
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

	/* Create the output for the list of microstates */
	if (nlhs >= 1)
	{
		plhs[0] = mxCreateDoubleMatrix(Ne, Nu, mxREAL);
		memcpy(mxGetPr(plhs[0]), microstates, Ne*Nu*sizeof(double));
	}

	/* Create the output for the classification of the signal */
	/* <!> indexes have to be converted to Matlab's convention ! */
	if (nlhs >= 2)
	{
		iSM = signalMicrostate;
		for (t = 0 ; t < Nt ; t++)
			(*iSM++) += 1;
		plhs[1] = mxCreateNumericMatrix(Nt, 1, mxUINT8_CLASS, mxREAL);
		memcpy(mxGetPr(plhs[1]), signalMicrostate, Nt*sizeof(char));
	}

	/* Create the output for the magnitude */
	if (nlhs >= 3)
	{
		plhs[2] = mxCreateDoubleMatrix(Nt, 1, mxREAL);
		memcpy(mxGetPr(plhs[2]), magnitude, Nt*sizeof(double));
	}

	/* Create the output for the measure R2 */
	if (nlhs >= 4)
		plhs[3] = mxCreateDoubleScalar(R2);

	/* Create the output for Sigmau2 */
	if (nlhs >= 5)
		plhs[4] = mxCreateDoubleScalar(sigmau2);

	/* Relesase memory */
	free(signal);
	free(microstates);
	free(magnitude);
	free(signalMicrostate);
}
