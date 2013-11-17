#ifndef _OMP_H_
#define _OMP_H_

#include "matrix.h"
#include "f2c.h"
#include "clapack.h"

using namespace std;

/****************************************************************************
Usage:
	Solving the sparse representation problem:
		min ||A*x - b||_2 + lambda * ||x||_1

Input:
	sols			  solution of OMP.
	y		          vector of length n.
	Dictionary        	  An explicit nxN matrix, as a dictionary.
        nSolver           	  length of solution vector. 
	nMaxIters		  maximum number of iterations to perform. 
        lambdaStop		  If specified, the algorithm stops when the last coefficient 
			   	  entered has residual correlation <= lambdaStop. 
	OptTolerance      	  Error tolerance, default 1e-5
*****************************************************************************/

bool solveOMP(vector<double> &sols,
			  vector<double> y,
			  Matrix Dictionary,
			  int nSolver,
			  int nMaxIters,
			  int lambdaStop = 0,
			  double OptTolerance = 1e-5);

Matrix linsolve(Matrix A, Matrix b);

bool updateChol(Matrix &_R,
				int _n,
				int _N,
				Matrix _dictionary,
				vector<int> _activeset,
				int _newindex);

#endif 
