#include "omp.h"

bool solveOMP(vector<double> &sols,
	vector<double> y,
	Matrix Dictionary,
	int nSolver,
	int nMaxIters,
	int lambdaStop,
	double OptTolerance)
{
	bool flag = false;	//for returning
	vector<vector<double>> temp;
	int n = y.size();

	Matrix _y(1,y.size());
	_y.setdata(vector<vector<double>>(1,y));
	_y = _y.trans();
	double normy = _y.norm();
	double resnorm = normy;

	Matrix x = zeros(nSolver,1);
	int k = 1;
	bool done = false;
	Matrix res(1,n);
	temp.assign(1,y);
	res.setdata(temp);
	res = res.trans();

	vector<int> Activeset;
	while(!done){
		/****************************************************************************
		Matlab Code:
		corr = Dic'*res;
		[maxcorr i] = max(abs(corr));
		****************************************************************************/
		Matrix corr = Dictionary.trans()*res;
		int newindex;
		double maxcorr = corr.m_abs().max_vec(newindex);
		Matrix R_I(1,1);

		//Update Cholesky factorization
		bool flag = updateChol(R_I,n,nSolver,Dictionary,Activeset,newindex);
		Activeset.push_back(newindex);
		
		//Solve for the least squares update.
		/****************************************************************************
		Matlab Code:
		dx = zeros(N,1);
		z = linsolve(R_I,corr(activeSet),opts_tr);
		dx(activeSet) = linsolve(R_I,z,opts);
		x(activeSet) = x(activeSet) + dx(activeSet);
		****************************************************************************/
		Matrix dx = zeros(nSolver,1);
		Matrix z = linsolve(R_I,corr.select(Activeset,1));
		Matrix dxactive = linsolve(R_I,z);
		for(int i=0;i<dxactive.Getcols();i++)
			dx(Activeset[i],1) = dxactive(i,1);
		for(vector<int>::size_type i=0;i<Activeset.size();i++)
			x(Activeset[i],1) = x(Activeset[i],1) + dx(Activeset[i],1);
		
		//Compute new residual
		/****************************************************************************
		Matlab Code:
		res = y - Dic(:,activeSet) * x(activeSet);
		****************************************************************************/
		Matrix residue = _y - Dictionary.select(Activeset,2) * x.select(Activeset,1);
		resnorm = residue.norm();
		if((resnorm<OptTolerance*normy||resnorm==OptTolerance*normy)||((lambdaStop>0)&&(maxcorr<lambdaStop||maxcorr==lambdaStop)))
			done = true;
		k++;
		if(k>nMaxIters)
			done = true;
		if(done){
			for(int i=0;i<x.Getrows();i++)
				sols[i] = x(i,1);
		}
	}
	return flag;
}

/****************************************************************************
Updates the Cholesky factor R of the matrix Dic(:,activeSet)' Dic(:,ActiveSet)
by adding Dic(:,newIndex) If the candidate column is in the span of the active 
set, R is not updated, and flag is set to 1.
****************************************************************************/
bool updateChol(Matrix &_R,
	int _n,
	int _N,
	Matrix _dictionary,
	vector<int> _activeset,
	int _newindex){
		bool flag = false;
		Matrix newVec = _dictionary.select(vector<int>(1,_newindex),2);
		if(!_activeset.size()){
			assert(_R.Getcols() == 1 && _R.Getrows() == 1);
			double temp = newVec.m_power(2).Vecsum(1);
			_R(0,0) = sqrt(temp);
		}
		//linear equation solver
		Matrix p = linsolve(_R.trans(),_dictionary.select(_activeset,2).trans()*(_dictionary.select(vector<int>(1,_newindex),2)));
		double q = newVec.m_power(2).Vecsum(1) + p.m_power(2).Vecsum(1);
		if(q < 1e-5 && q == 1e-5)
			flag = true;
		else{
			//pooling four matrix
			/****************************************************************************
			Matlab Code:
			_R = [_R p; zeros(1, size(_R,2)) sqrt(q)];
			****************************************************************************/
			Matrix res(_R.Getrows(),_R.Getcols());
			for(int i=0;i<_R.Getrows();i++)
				for(int j=0;j<_R.Getcols();j++)
					res(i,j) = _R(i,j);
			for(int i=0;i<p.Getcols();i++)
				res(i,_R.Getcols()) = p(i,1);
			res(res.Getrows(),res.Getcols()) = sqrt(q);
			_R = res.copy();
		}
		return flag;
}

/****************************************************************************
This function arims at the linear quation solver and the CLAPACK is in need.
I adopt the function of dgesv_, which deals with case of A being square.
****************************************************************************/
Matrix linsolve(Matrix A, Matrix b){
	//A*X = b solver
	assert(A.Getrows() == b.Getrows());
	//Solving the case of A being square!
	assert(A.Getcols() == A.Getrows());

	Matrix res(A.Getcols(),b.Getrows());
	integer i, j, N, nrhs, lda, ldb, info;
	N = A.Getcols();
	lda = N;
	ldb = N;
	nrhs = b.Getcols();

	integer *ipiv = new integer[N]();
	doublereal *_A = new doublereal[N*N]();
	doublereal *_B = new doublereal[N*nrhs]();

	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			_A[i*N+j] = A(j,i);
	for(i=0;i<N;i++)			
		for(j=0;j<nrhs;j++)
			_B[i*N+j] = b(j,i);

	dgesv_(&N,&nrhs,_A,&lda,ipiv,_B,&ldb,&info);

	if(!info){
		cout<<"Linear Equation Solving FAILED!"<<endl;
		exit(-1);
	}

	for(i=0;i<N;i++)		
		for(j=0;j<nrhs;j++)
			res(j,i) = _B[i*N+j];

	return res;
}

