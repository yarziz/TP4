#ifndef _LAPLACIAN_CPP

#include "Laplacian.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
Laplacian::Laplacian(Function* function, DataFile* data_file) :
  _fct(function), _df(data_file)
{
}

void Laplacian::InitialCondition(VectorXd & sol){
  int Nx=_df->Get_Nx();
  int Ny=_df->Get_Ny();
  double hx=_df->Get_hx();
  double hy=_df->Get_hy();
  int Nx_Ny=Nx*Ny;
  VectorXd sol(Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      sol[j*Nx+i]=_fct->_fct->InitialCondition(_df->Get_xmin()+(i+1)*hx,_df->Get_ymin()+(j+1)*hy);
    }
  }
}

void Laplacian::BuildLaplacianMatrix(){
  int Nx=_df->Get_Nx();
  int Ny=_df->Get_Ny();
  double hx=_df->Get_hx();
  double hy=_df->Get_hy();
  int Nx_Ny=Nx*Ny;
  // SparseMatrix<double> H(Nx_Ny,Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      if(i<Nx-1){
	_H.insert(j*Nx+i,j*Nx+i+1,-1/(hx*hx));
      }
      if(i>0){
	_H.insert(j*Nx+i,j*Nx+i-1,-1/(hx*hx));
      }
      _H.insert(j*Nx+i,j*Nx+i+,2/(hx*hx)+2/(hy*hy));
      if(j>0){
	_H.insert(j*Nx+i,(j-1)*Nx+i,-1/(hy*hy));
      }
      if(j<Ny-1){
	_H.insert(j*Nx+i,(j+1)*Nx+i,-1/(hy*hy));
      }
    }
  }
}

void Laplacian::BuildSourceTerm(VectorXd & f,double t){
  int Nx=_df->Get_Nx();
  int Ny=_df->Get_Ny();
  double hx=_df->Get_hx();
  double hy=_df->Get_hy();
  int Nx_Ny=Nx*Ny;
  VectorXd sol(Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      f[j*Nx+i]=_fct->_fct->SourceFunction(_df->Get_xmin()+(i+1)*hx,_df->Get_ymin()+(j+1)*hy,t);
    }
  }
}


VectorXd Laplacian::ExactSolution(double t){
  int Nx=_df->Get_Nx();
  int Ny=_df->Get_Ny();
  double hx=_df->Get_hx();
  double hy=_df->Get_hy();
  int Nx_Ny=Nx*Ny;
  VectorXd sol(Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      sol[j*Nx+i]=_fct->_fct->ExactSolution(_df->Get_xmin()+(i+1)*hx,_df->Get_ymin()+(j+1)*hy,t);
    }
  }
  return sol;
}


void Laplacian::BuildF(const VectorXd & sol,const VectorXd& f){
  _p=(-Get_sigma()*(_H*sol)+f);
}




#define _LAPLACIAN_CPP
#endif
