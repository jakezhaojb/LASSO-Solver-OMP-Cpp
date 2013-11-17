//Designed by Junbo ZHAO
//2013.11.15

#include "omp.h"
#include "matrix.h"

void main(){
	Matrix Dic;
	Dic.readfile("D.txt");
	ifstream finy("y.txt");
	assert(finy);
	vector<double> y;
	double temp = 0;
	while(finy.peek()!=EOF){
		finy>>temp;
		y.push_back(temp);
	}
	finy.close(); 
	vector<double> sols;
	//OMP running with a 150 dim sparse solution and 500 iteraions requested.
	if(solveOMP(sols,y,Dic,150,500))
		cout<<"OMP for LASSO accomplished"<<endl;
	ofstream fout("solution.txt");
	fout.setf(ios::left);
	for(vector<double>::size_type i=0;i!=sols.size();i++){
		fout.width(15);
		fout<<sols[i]<<" "<<endl;
	}
}
