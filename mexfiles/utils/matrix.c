/*
	matrix.c

	Copyright (C) 2012 Frederic Wilhelm <frederic.wilhelm at epfl.ch>
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

/*************************************************************************
 *
 *		INCLUDES
 *
 *************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tred2.h"
#include "tqli.h"
#include "matrix.h"

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

/* \brief Transpose a matrix
 *
 * \params  matrix
 *		Matrix pointer (content is modified !)
 *	\params	m
 *		Number of lines (Matlab convention, i.e. matrix[n][m])
 *	\params	n
 *		Number of columns (Matlab convention, i.e. matrix[n][m]);
 *	\return
 *		Pointer to the matrix (i.e. param matrix itself)
 */
double *transpose(double *matrix, int m, int n)
{
	double *buffer = malloc(m*n*sizeof(double));
	double *iB = buffer;
	int i,j;

	/* Copy matrix into buffer */
	memcpy(buffer, matrix, m*n*sizeof(double));

	/* Transpose */
	for (j=0 ; j<n ; j++)
		for (i=0 ; i<m ; i++)
		{
			*(matrix + i*n + j) = *iB++;
		}

	/* Free buffer */
	free(buffer);

	/* Return pointer to matrix */
	return matrix;
}

/* \brief Finds the largest eigenvector of a symmetric matrix
 *
 * Uses Householder reduction (code taken from Numerical Recipes in C)
 *
 * \params  matrix
 *		Matrix pointer (content is modified !)
 *	\params	n
 *		Matrix dimension
 *	\params	e
 *		Largest eigenvector
 */
void largestEigenvector(double *matrix, int n, double *eigenvector)
{
	/* <!> *$€%+§&@ convention of Numerical Recipes to begin index by 1 */

	/* Copy matrix into n*n table */
	int i;
	double **a, *iM;
	double d[n+1]; /* diagonal elements */
	double e[n+1]; /* off-diagonal elements */
	double max; /* maximum eigenvalue */
	int iMax; /* index of maximum eigenvector */

	/* Initialize a and z */
	a = malloc((n+1)*sizeof(double *));
	iM = matrix;
	for (i = 1 ; i <= n ; i++)
	{
		a[i] = malloc((n+1)*sizeof(double));
		memcpy(a[i]+1, iM, n*sizeof(double));
		iM += n;
	}

	/* Householder reduction */
	tred2(a, n, d, e);

	/* QL algorithm */
	tqli(d, e, n, a);

	/* Search eigenvector of maximum eigenvalue */
	max = 0;
	iMax = 1;
	for (i = 1 ; i <= n ; i++)
	{
		if (d[i] > max)
		{
			max = d[i];
			iMax = i;
		}
	}

	/* Return maximum eigenvector */
	max = 0;
	for (i = 0 ; i < n ; i++)
	{
		eigenvector[i] = *(a[i+1] + iMax);
		max += SQR(eigenvector[i]);
	}
	max = sqrt(max);
	for (i = 0 ; i < n ; i++)
		eigenvector[i] /= sqrt(max);

	/* release space */
	for (i = 1 ; i <= n ; i++)
	{
		free(a[i]);
	}
	free(a);
}
