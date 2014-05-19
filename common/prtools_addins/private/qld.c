	/*
 * QLD.C 
 *
 * A Matlab MEX interface to Prof. K. Schittkowski's (Univ. of Bayreuth,
 * Germany) and Prof. M.J.D. Powell's (Univ. of Cambridge, UK) FORTRAN QLD 
 * routine, kindly provided by Prof. A.L. Tits (Univ. of Maryland, US).
 *
 * This routine can be called as:
 *
 * [x, u] = QLD (H, f, A, b, lb, ub, x0, n)
 *
 * where 
 *
 * x	= vector that minimizes: 0.5*x' H x + f' x,
 *
 *      with constraints A(i)x + b(i)  = 0, 	i = 1, ..., n
 *                       A(i)x + b(i) >= 0,   i > n (????)
 *
 *                   and lb < x < ub 
 *        
 * u  = set of Lagrange multipliers of the constraints, lower bounds
 *      and upper bounds (in that order) in point x
 *
 * The function will return x = [] and print an apropriate message in
 * case an error occurs.
 *
 * x0 can be used as an initial guess (not implemented yet).
 *
 * (C) 1996 David Tax & Dick de Ridder
 * Faculty of Applied Physics, Delft University of Technology
 * P.O. Box 5046, 2600 GA Delft, The Netherlands 
 *
 * $Id: qld.c,v 1.1.1.1 2006/03/07 15:20:39 davidt Exp $
 */
 
#include <stdio.h>
#include <string.h>
#include <mex.h>

#undef DEBUG

int ql0001_(int *m,int *me,int *mmax,int *n,int *nmax,int *mnn,
            double *c,double *d,double *a,double *b,double *xl,
            double *xu,double *x,double *u,int *iout,int *ifail,
            int *iprint,double *war,int *lwar,int *iwar,int *liwar,
            double *eps1);
                                                
#define	VAR_H		(0)
#define	VAR_f		(1)
#define VAR_A		(2)
#define VAR_b		(3)
#define VAR_vlb	(4)
#define VAR_vub	(5)
#define VAR_xo	(6)
#define VAR_n		(7)

#define Matrix mxArray

void printMyMatrix (char *name, Matrix *m)
{
	int			 i, j;
	double	*p;

	p = mxGetPr (m);
	mexPrintf ("%s = \n", name);
	for (i = 1; i <= mxGetN (m); i++)
	{
		for (j = 1; j <= mxGetM (m); j++)
		{
			mexPrintf ("%lf ", *p);
			p++;
		}
		mexPrintf ("\n");
	}
}

void mexFunction (int nlhs, Matrix *plhs[],
          				int nrhs, const Matrix *prhs[])
{
	Matrix	*x, 
					*u,
					*error;
					
	int			*no_constraints,
					*no_equality_constraints,
					*no_variables,
					*no_mnn,
					*output_fd,
					*output_fail,
					*output_print,
					*no_double_working_array,
					*no_int_working_array,
					*int_working_array;
					
	double	*H,
					*f,
					*A,
					*b,
					*vlb,
					*vub,
					*epsilon,
					*double_working_array;
						
	/* We need 8 arguments, and return 1 or 2. */
	
	if ((nrhs != 8) || (nlhs < 1) || (nlhs > 2))
	{
		mexPrintf ("usage: [x, u] = QLD (H,f,A,b,lb,ub,x0,n)\n");
		return;
	}
	
	no_constraints = (int *) malloc (sizeof (int));
		*no_constraints = mxGetM (prhs[VAR_A]);

	no_equality_constraints = (int *) malloc (sizeof (int));
		*no_equality_constraints = mxGetScalar (prhs[VAR_n]);
		
	no_variables = (int *) malloc (sizeof (int));
		*no_variables = mxGetN (prhs[VAR_A]);
		
	no_mnn = (int *) malloc (sizeof (int));
		*no_mnn = *no_constraints + 2 * *no_variables;

	H 		= (double *) mxGetPr (prhs[VAR_H]);
	f 		= (double *) mxGetPr (prhs[VAR_f]);
	A 		= (double *) mxGetPr (prhs[VAR_A]);
	b 		= (double *) mxGetPr (prhs[VAR_b]);
	vlb 	= (double *) mxGetPr (prhs[VAR_vlb]);
	vub 	= (double *) mxGetPr (prhs[VAR_vub]);
	
	x			= mxCreateDoubleMatrix (*no_variables, 1, mxREAL);
	u			= mxCreateDoubleMatrix (*no_mnn,       1, mxREAL);
	error = mxCreateDoubleMatrix (0,             0, mxREAL);
	
	no_double_working_array = (int *) malloc (sizeof (int));
		*no_double_working_array = (int) ((3 * *no_variables * *no_variables) / 2 + 10 * *no_variables + 2 * *no_constraints + 1);

	no_int_working_array = (int *) malloc (sizeof (int));	
		*no_int_working_array = *no_variables;
		
	double_working_array = (double *) malloc (*no_double_working_array * sizeof (double)); 

	int_working_array = (int *) malloc (*no_int_working_array * sizeof (int));
		int_working_array[0] = 1; int_working_array[1] = 1;
	
	output_fd = (int *) malloc (sizeof (int));
		*output_fd = 1;
	output_fail = (int *) malloc (sizeof (int));
	output_print = (int *) malloc (sizeof (int));
		*output_print = 1;	
	
	epsilon = (double *) malloc (sizeof (double));
		*epsilon = 1.0e-20;

#ifdef DEBUG
 	mexPrintf ("No. constraints         : %d\n", *no_constraints);
 	mexPrintf ("No. equality constraints: %d\n", *no_equality_constraints);
 	mexPrintf ("No. variables           : %d\n", *no_variables);

	printMyMatrix ("H", prhs[VAR_H]);
	printMyMatrix ("f", prhs[VAR_f]);
	printMyMatrix ("A", prhs[VAR_A]);
	printMyMatrix ("b", prhs[VAR_b]);
	printMyMatrix ("vlb", prhs[VAR_vlb]);
	printMyMatrix ("vub", prhs[VAR_vub]);
#endif

	ql0001_ (no_constraints,
					 no_equality_constraints,
					 no_constraints,
					 no_variables,
					 no_variables,
					 no_mnn,
					 H, f, A, b, vlb, vub, 
					 mxGetPr (x),
					 mxGetPr (u),
					 output_fd,
					 output_fail,
					 output_print,
					 double_working_array,
					 no_double_working_array,
					 int_working_array,
					 no_int_working_array,
					 epsilon);

	if (*output_fail > 0)
	{
		plhs[0] = error;
		
		if (nlhs > 1)
		{
			plhs[1] = error;
		}
		
/*		mxFree (x);
		mxFree (u);*/
    mxDestroyArray(x);
    mxDestroyArray(u);
			
   	if (*output_fail == 1)
   	{
   		mexPrintf ("error: too many iterations\n");
   	}
   	else if (*output_fail == 2)
   	{
   		mexPrintf ("error: insufficient accuracy to satisfy convergence "
   			         		"criterion\n");
   	}
   	else if (*output_fail == 5)
   	{
   		mexPrintf ("error: internal - working array too short\n");
   	}
   	else if (*output_fail > 10)
   	{
   		mexPrintf ("error: inconsistent constraints\n");
   	}
	}		
	else
	{
		plhs[0] = x;
		
		if (nlhs == 2)
		{
			plhs[1] = u;
		}
		else
		{
/*			mxFree (u);*/
      mxDestroyArray(u);
		}
			
#ifdef DEBUG
		printMyMatrix (x);
#endif	
	}								 

	free (epsilon);
	free (output_print);
	free (output_fail);
	free (output_fd);
	free (int_working_array);
	free (double_working_array);
	free (no_int_working_array);
	free (no_double_working_array);
	free (no_mnn);
	free (no_variables);
	free (no_equality_constraints);
	free (no_constraints);
}
