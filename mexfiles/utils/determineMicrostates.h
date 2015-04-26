/*
	matrix.h

	Copyright (C) 2012 Frederic Wilhelm <frederic.wilhelm at epfl.ch>
*/

/*************************************************************************
 *
 *		FUNCTIONS
 *
 *************************************************************************/

/* Determine the microstates of a time samples */
void determineMicrostates(int Nu,
		int Ne,
		int Nt,
		double *signal, /* (Ne * Nt) */
		double *microstates, /* (Ne * Nu) */
		char *signalMicrostate, /* Nt */
		double* magnitude); /* Nt */

/* Determines the responsibilities of a time samples w.r.t microstates */
void computeResponsibilities(int Nu,
		int Ns,
		int Nt,
		double *signal, /* (Ne * Nt) */
		double *microstates,  /* (Ne * Nu) */
		double *sigma2, /* (Nu) */
		double *mix, /* (Nu) */
		double *gamma); /* (Ne * Nt) */
