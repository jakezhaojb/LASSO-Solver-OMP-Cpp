#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <algorithm>

using namespace std;

class Matrix{
private:
	int rows;
	int cols;
	vector<vector<double>> data;
public:
	vector<double> eigen;
	vector<vector<double>> eigenvec;
	Matrix(){
		rows = 0;
		cols = 0;
	}
	Matrix(int r,int c){
		rows = r;
		cols = c;
		data.assign(rows,vector<double>(cols,0));
	}
	void readfile(string fname);
	void writefile(string fname);
	Matrix operator*(Matrix &mul) const;
	Matrix operator+(Matrix &ad) const;
	Matrix operator-(Matrix &su) const;
	Matrix operator/(const double &div) const;
	Matrix trans();
	Matrix m_sqrt();
	vector<double> mean(int option);
	Matrix sort(vector<int>& id,int option); 
	Matrix trunc(int num, int option);
	Matrix invert();
	Matrix copy(){
		Matrix res(rows,cols);
		res.setdata(data);
		return res;
	}
	double& operator()(int x,int y){
		return data[x][y];
	}
	const double& operator()(int x,int y) const{
		return data[x][y];
	}

	Matrix select(vector<int> num,int option);

	int Getrows() const{
		return this->rows;
	}

	int Getcols() const{
		return this->cols;
	}

	vector<vector<double>> Getdata() const{
		return this->data;
	}

	void setdata(vector<vector<double>> d){
		if(data.size()!=d.size()){
			cout<<"Matrix dimensions must agree."<<endl;
			exit(-1);
		}
		data.assign(d.begin(),d.end());
	}

	bool symm() const{
		if(rows!=cols)
			return false;
		else{
			for(int i=0;i<rows;i++)
				for(int j=0;j<cols;j++)
					if(abs(data[i][j]-data[j][i])>1e-3)
						return false;
			return true;
		}
	}
	double max_vec(int &l);
	Matrix m_abs();
	Matrix m_power(int p);
	
	double Vecsum(int option){
		double s;
		if(option == 1){//rows
			assert(cols==1);
			s = 0;
			for(int i=0;i<rows;i++)
				s = s + data[i][1];
		}
		else if(option == 2)//columns
			trans().Vecsum(1);
		else{
			cerr<<"WRONG Input arguments"<<endl;
			exit(-1);
		}
		return s;
	}
	double norm();
};
inline Matrix zeros(int r,int c){
	Matrix res(r,c);
	return res;
}
inline Matrix ones(int r,int c){
	Matrix res(r,c);
	vector<vector<double>> data_t(r,vector<double>(c,1));
	res.setdata(data_t);
	return res;
}
inline Matrix diag(vector<double> d){
	int n = d.size();
	Matrix res(n,n);
	vector<vector<double>> data_t(n,vector<double>(n,0));
	for(int i=0;i<n;i++)
		data_t[i][i] = d[i];
	res.setdata(data_t);
	return res;
}

#endif
