/*
	gaussian.c

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

#include <math.h>
#include "gaussian.h"

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

/* \brief Gaussian function
 *
 * \params  V
 *		Main parameter (signal)
 *	\params	u
 *		Microstate
 *	\params	sigma2
 *		Variance
 * \params	Ns
 *		Dimension (number of channels)
 */
double gaussian(double* V, double* u, double sigma2, int Ns)
{
	int i;
	double *iV, *iU;
	double norm, scal;

	/* Distance of the signal to the microstate */
	iV = V;
	iU = u;
	norm = scal = 0;
	for (i = 0 ; i < Ns ; i++)
	{
		norm += SQR(*iV);
		scal += (*iV++) * (*iU++);
	}
	norm = norm - scal*scal;

	/* Calculate value */
	return exp(-norm/(2*sigma2)) / pow(sqrt(2*M_PI*sigma2), Ns-2);
}
