/*
	matrix.h

	Copyright (C) 2012 Frederic Wilhelm <frederic.wilhelm at epfl.ch>
*/

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* Transpose a matrix */
double *transpose(double *matrix, int m, int n);

/* Finds the largest eigenvector of a symmetric matrix */
void largestEigenvector(double *matrix, int n, double *eigenvector);
