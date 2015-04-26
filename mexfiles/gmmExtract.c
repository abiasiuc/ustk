/*
	gmmExtract.c

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
#include "utils/gaussian.h"
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

/* \brief MEX routine for microstate extraction using GMM
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
 *  \params prhs[3]
 *      (optional) Impose same variance for all microstate (default 0)
 *	\params	prhs[4]
 *		(optional) Random seed
 *	\params	prhs[5]
 *		(optional) Initial microstates
 *	\params	prhs[6]
 *	   (optional) Initial variances
 *	\params	prhs[7]
 *		(optional) Initial mixing coefficients (1 * Nu)
 * 
 *	\return	plhs[0]
 *		List of states (number of electrodes * number of states)
 *	\return	plhs[1]
 *		Responsibilities (time * Nu)
 *	\return	plhs[2]
 *		Variances of each microstate (1 * Nu)
 *	\return	plhs[3]
 *		Fit measure R^2
 *	\return	plhs[4]
 *		Estimated noise variance Sigmau2
 *	\return	plhs[5]
 *		Mixing coefficients (1 * Nu)
 */ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int t, i, k, l;
	int Nu, Ns, Nt;
	double norm, scal;
	double *signal, *microstates, *iS, *iS2, *iM;
	double *sigma2, *iSig, sigmau2; /* Variances */
	double *mix, *iMix ; /* Mixing coefficients */
	double *gamma, *iG; /* Responsibilities */
	double *matS, *iMatS;
	double Nk;
	double logLik, logLik0=0, deviation;
	double epsilon; /* Radius of convergence */
	bool stop=false; /* stop the loop */
   bool cstSigma2 = false; /* Use a constant sigma */
    
	if (nrhs < 3)
		mexErrMsgTxt("Must have at least three input argument");		

	/* Retrieve the seed */
	if (nrhs >= 5)
		srand((int) *mxGetPr(prhs[4]));
	else
		srand(time(0));

	/* Retrieve the signal data AND TRANSPOSE IT (for later convenience) */
	Nt = mxGetM(prhs[0]);
	Ns = mxGetN(prhs[0]);
	signal = malloc(Nt * Ns * sizeof(double));
	memcpy(signal, mxGetPr(prhs[0]), Ns*Nt*sizeof(double));
	transpose(signal, Nt, Ns);

	/* Retrieve the number of microstates */
	Nu = (int) *mxGetPr(prhs[1]);

	/* Retrieve the radius of convergence */
	epsilon = *mxGetPr(prhs[2]);
    
    /* Retrieve the option to use a unique variance */
    if (nrhs >= 4 && fabs(*mxGetPr(prhs[3])) > 0.0001)
        cstSigma2 = true;
	
	/* Initialize the microstates */
	microstates = malloc(Nu * Ns * sizeof(double));
	if (nrhs >= 6)
		memcpy(microstates, mxGetPr(prhs[5]), Nu*Ns*sizeof(double));
	else
	{
		/* Initialize iterator */
		iM = microstates;
		for (i = 0 ; i < Nu ; i++)
		{
			/* Generate one random microstate */
			scal = 0; /* Average (for transformation to ARD) */
			norm = 0; /* Norm (for normalization) */
			for (k = 0 ; k < Ns ; k++)
			{
				*iM = rand() % 1000;
				scal += *iM;
				iM++;
			}
			/* Transform to average reference data */
			iM -= Ns;
			scal /= Ns;
			for (k = 0 ; k < Ns ; k++)
			{
				*iM -= scal;
				norm += SQR(*iM);
				iM++;
			}
			/* Normalize the microstate */
			iM -= Ns;
			norm = sqrt(norm);
			for (k = 0 ; k < Ns ; k++)
			{
				*iM++ /= norm;
			}
		}
	}

	/* Initialize the variances */
	sigma2 = malloc(Nu * sizeof(double));
	if (nrhs >= 7)
		memcpy(sigma2, mxGetPr(prhs[6]), Nu*sizeof(double));
	else
	{
		iSig = sigma2;
		for (i = 0 ; i < Nu ; i++) (*iSig++) = 1.0;
	}

	/* Initialize the mixing coefficients */
	mix = malloc(Nu * sizeof(double));
	if (nrhs >= 8)
		memcpy(mix, mxGetPr(prhs[7]), Nu*sizeof(double));
	else
	{
		iMix = mix;
		for (i = 0 ; i < Nu ; i++) (*iMix++) = 1.0/Nu;
	}

	/* Initialize the responsibilities */
	gamma = malloc(Nt*Nu*sizeof(double));

	/* Initialize the S matrices */
	matS = malloc(Ns * Ns * sizeof(double));

	/* Main loop */
	int count = 0;
	while (!stop)
	{
		/* E STEP */
		/* Compute the responsibilities */
		computeResponsibilities(Nu, Ns, Nt, signal, microstates, sigma2,
				mix, gamma);

		/* M STEP */
		iM = microstates;
		iSig = sigma2;
		iMix = mix;
		sigmau2 = 0;
		for (i = 0 ; i < Nu ; i++)
		{
			/* Compute matrix S and Nk */
			Nk = 0;
			iG = gamma + i;
			iS = signal;
			memset(matS, 0, Ns*Ns*sizeof(double));
			for (t = 0 ; t < Nt ; t++)
			{
				Nk += *iG;
				iMatS = matS;
				iS2 = iS;
				for (k = 0 ; k < Ns ; k++)
				{
					for (l = 0 ; l < Ns ; l++)
					{
						(*iMatS++) += (*iG) * (*iS) * (*iS2++);
					}
					iS++;
					iS2 -= Ns;
				}

				/* Next signal */
				iG += Nu;
			}

			/* Extract microstate */
			largestEigenvector(matS, Ns, iM);

			/* Extract variance */
			*iSig = 0;
			iS = signal;
			iG = gamma + i;
			for (t = 0 ; t < Nt ; t++)
			{
				scal = 0;
				for (k = 0 ; k < Ns ; k++)
				{
					*iSig += (*iG) * (*iS) * (*iS);
					scal += (*iS++) * (*iM++);
				}
				*iSig -= (*iG) * SQR(scal);

				/* Next signal */
				iG += Nu;
				iM -= Ns;
			}
			*iSig = fmax(*iSig, 1e-5*Nk*(Ns-2)); /* To avoid kaboom */
			sigmau2 += (*iSig) / (Nt*(Ns-1));
			*iSig++ /= (Nk*(Ns-2));

			/* Mixing coefficient */
			*iMix++ = Nk/Nt;

			/* Next microstate */
			iM += Ns;
		}

		/* Same sigma for all map (if option is set to true) */
        if (cstSigma2)
            for (i = 0 ; i < Nu ; i++)
            {
                sigma2[i] = sigmau2;
            }

		/* Log likelihood */
		logLik = 0;
		iS = signal;
		for (t = 0 ; t < Nt ; t++)
		{
			/* Compute the log for each gaussian */
			scal = 0; /* temp value */
			iMix = mix;
			iM = microstates;
			iSig = sigma2;
			for (i = 0 ; i < Nu ; i++)
			{
				scal = scal + (*iMix++) * gaussian(iS, iM, *iSig++, Ns);
				iM += Ns;
			}

			/* Updata log likelihood */
			logLik += log(scal);

			/* Next signal */
			iS += Ns;
		}

		/* Compute deviation */
		deviation = fabs(logLik - logLik0)/(fabs(logLik)+fabs(logLik0));
		if (deviation < epsilon)
		{
			stop = true;
		}
		else
		{
			logLik0 = logLik;
		}
	}

	/* Create the output for the list of microstates */
	if (nlhs >= 1)
	{
		plhs[0] = mxCreateDoubleMatrix(Ns, Nu, mxREAL);
		memcpy(mxGetPr(plhs[0]), microstates, Ns*Nu*sizeof(double));
	}

	/* Create the output for the responsibilities */
	if (nlhs >= 2)
	{
		plhs[1] = mxCreateDoubleMatrix(Nt, Nu, mxREAL);
		computeResponsibilities(Nu, Ns, Nt, signal, microstates, sigma2,
				mix, mxGetPr(plhs[1])); /* Update gamma (bcs of the last iteration) */
		transpose(mxGetPr(plhs[1]), Nu, Nt);
	}

	/* Create the output for the variances */
	if (nlhs >= 3)
	{
		plhs[2] = mxCreateDoubleMatrix(Nu, 1, mxREAL);
		memcpy(mxGetPr(plhs[2]), sigma2, Nu*sizeof(double));
	}

	/* Create the output for R^2 */
	if (nlhs >= 4)
	{
		/* Compute sigmad2 (sum of (Vt)^2) */
		iS = signal;
		norm = 0;
		for (k = 0 ; k < Ns*Nt ; k++)
		{
			norm += SQR(*iS);
			iS++;
		}
		norm = norm / (Nt * (Ns - 1)); /* sigmad2 */

		/* return R^2 */
		plhs[3] = mxCreateDoubleScalar(1 - sigmau2/norm);
	}

	/* Create the output for sigmau2 */
	if (nlhs >= 5)
	{
		plhs[4] = mxCreateDoubleScalar(sigmau2);
	}

	/* Create the output for the mixing coefficients) */
	if (nlhs >= 6)
	{
		plhs[5] = mxCreateDoubleMatrix(Nu, 1, mxREAL);
		memcpy(mxGetPr(plhs[5]), mix, Nu*sizeof(double));
	}

	/* Relesase memory */
	free(signal);
	free(microstates);
	free(gamma);
	free(sigma2);
	free(mix);
	free(matS);
}

