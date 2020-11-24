#ifndef _LAPLACIAN_H

#include <string>
#include "Dense"
#include "Sparse"
#include "Function.h"

class Laplacian
{
private:
	Function* _fct;
	DataFile* _df;
       	Eigen::VectorXd _p;
	Eigen::SparseMatrix<double> _H;

public:
	// Constructeur
        Laplacian(Function* function, DataFile* data_file);
	void InitialCondition(Eigen::VectorXd & sol);
	void BuildLaplacianMatrix();
	void BuildSourceTerm(Eigen::VectorXd & f,double t);
	void BuildF(Eigen::VectorXd sol,Eigen::VectorXd f);
	Eigen::VectorXd ExactSolution(double t);
	const Eigen::VectorXd &  Get_p()  const;
	

};

#define _LAPLACIAN_H
#endif
