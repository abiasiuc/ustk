/*
	determineMicrostates.c

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
#include "gaussian.h"
#include "determineMicrostates.h"

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* \brief Determine the microstates of time samples
 *
 * "Optional" parameters means that they are pointers which can be replaced
 * by NULL.
 *
 * \params  Nu
 *		Number of microstates
 *	\params	Ns
 *		Number of electrodes
 *	\params	Nt
 *		Number of time samples
 * \params	signal
 *		Time samples (Ns * Nt)
 *	\params	microstates
 *		Microstates (Ns * Nu)
 *	\params	signalMicrostate
 *		Time serie of corresponding microstates (Nt)
 *	\parms	magnitude
 *		Time serie of magnitudes (Nt) (optional)
 */
void determineMicrostates(int Nu,
		int Ns,
		int Nt,
		double *signal, /* (Ns * Nt) */
		double *microstates, /* (Ns * Nu) */
		char *signalMicrostate, /* Nt */
		double* magnitude) /* Nt */
{
	int t, i, k;
	int sm;
	double magn, temp;
	double *iS = signal, *iM,
			 *iMagn = magnitude; /* iterators */
	char	 *iSM = signalMicrostate;

	/* Loop over each time */
	for (t=0 ; t<Nt ; t++)
	{
		/* Loop over the microstates */
		magn = 0;
		sm = 0;
		iM = microstates;
		for (i=0 ; i<Nu ; i++)
		{
			/* Compute scalar product between signal at t and microstate i */
			temp = 0;
			for (k=0 ; k<Ns ; k++)
			{
				temp += iS[k] * (*iM++); /* (microstates+i*Ns)[k] * (signal + k*Nt)[t]; */
			}
			if (temp * temp > magn * magn)
			{
				magn = temp;
				sm = i;
			}

		}

		/* Attribute microstate sm & magnitude magn to t */
		*iSM++ = (char) sm;
		if (magnitude != NULL)
			*iMagn++ = magn;

		/* Jump to next time */
		iS += Ns;
	}
}

/* \brief Determines the responsibilities of a time samples w.r.t microstates
 *
 * \params  Nu
 *		Number of microstates
 *	\params	Ns
 *		Number of electrodes
 *	\params	Nt
 *		Number of time samples
 * \params	signal
 *		Time samples (Ns * Nt)
 *	\params	microstates
 *		Microstates (Ns * Nu)
 *	\params sigma2
 *		Variances of the gaussians (Nu)
 *	\params mix
 *		Mixing coefficients of the gaussians (Nu)
 *	\params	gamma
 *		Responsibilities of the signals w.r.t the microstates (Ns * Nt)
 */
void computeResponsibilities(int Nu,
		int Ns,
		int Nt,
		double *signal, /* (Ns * Nt) */
		double *microstates,  /* (Ns * Nu) */
		double *sigma2, /* (Nu) */
		double *mix, /* (Nu) */
		double *gamma) /* (Nu * Nt) */
{
	double *iG = gamma, *iS = signal, *iMix, *iM, *iSig, norm;
	int t, i;

	/* Compute the responsibilities at each time */
	for (t = 0 ; t < Nt ; t++)
	{
		norm = 0;
		iMix = mix;
		iM = microstates;
		iSig = sigma2;
		for (i = 0 ; i < Nu ; i++)
		{
			*iG = (*iMix++) * gaussian(iS, iM, *iSig++, Ns);
			norm += *iG++;
			iM += Ns;
		}

		/* Normalization */
		iG -= Nu;
		for (i = 0 ; i < Nu ; i++)
		{
			*iG++ /= norm;
		}

		/* Next signal */
		iS += Ns;
	}
}
