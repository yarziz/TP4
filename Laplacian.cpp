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

Laplacian::InitialCondition(){
  int Nx=data_file->Get_Nx();
  int Ny=data_file->Get_Ny();
  double hx=data_file->Get_hx();
  double hy=data_file->Get_hy();
  int Nx_Ny=Nx*Ny;
  VectorXd sol(Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      sol[j*Nx+i]=_fct->InitialCondition(_df->Get_xmin()+(i+1)*hx,_df->Get_ymin()+(j+1)*hy);
    }
  }
}

Laplacian::BuildLaplacianMatrix(){
  int Nx=data_file->Get_Nx();
  int Ny=data_file->Get_Ny();
  double hx=data_file->Get_hx();
  double hy=data_file->Get_hy();
  int Nx_Ny=Nx*Ny;
  SparseMatrix<double> H(Nx_Ny,Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      if(i<Nx-1){
	H.insert(j*Nx+i,j*Nx+i+1,-1/(hx*hx));
      }
      if(i>0){
	H.insert(j*Nx+i,j*Nx+i-1,-1/(hx*hx));
      }
      H.insert(j*Nx+i,j*Nx+i+,2/(hx*hx)+2/(hy*hy));
      if(j>0){
	H.insert(j*Nx+i,(j-1)*Nx+i,-1/(hy*hy));
      }
      if(j<Ny-1){
	H.insert(j*Nx+i,(j+1)*Nx+i,-1/(hy*hy));
      }
    }
  }
}

Laplacian::BuildSourceTerm(double t){
  int Nx=data_file->Get_Nx();
  int Ny=data_file->Get_Ny();
  double hx=data_file->Get_hx();
  double hy=data_file->Get_hy();
  int Nx_Ny=Nx*Ny;
  VectorXd sol(Nx_Ny);
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      f[j*Nx+i]=_fct->SourceFunction(_df->Get_xmin()+(i+1)*hx,_df->Get_ymin()+(j+1)*hy,t);
    }
  }
}




#define _LAPLACIAN_CPP
#endif
