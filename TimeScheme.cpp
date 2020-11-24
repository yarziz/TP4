#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>
#include <fstream>
#include <math.h>


using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, Laplacian* lap) :
  _df(data_file), _lap(lap)
{
}

EulerScheme::EulerScheme(DataFile* data_file, Laplacian* lap):
TimeScheme(data_file, lap) 
{
  
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, Laplacian* lap):
TimeScheme(data_file, lap) 
{
  
}

// Destructeur
TimeScheme::~TimeScheme()
{
}

void TimeScheme::Initialise()
{
  _lap->InitialCondition(_sol);
  _lap->BuildSourceTerm(_f,_df->Get_t0());
  _t=_df->Get_t0();
  _lap->BuildLaplacianMatrix();
  _lap->BuildSourceTerm(_f,_df->Get_t0());
  _n=0;
}



void TimeScheme::SaveSolution(int n)
{
  string n_file = _df->Get_results() + to_string(n) + ".vtk";
  ofstream solution;
  solution.open(n_file, ios::out);
  int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
  double xmin(_df->Get_xmin()), ymin(_df->Get_ymin());
  double hx(_df->Get_hx()), hy(_df->Get_hy());

  solution << "# vtk DataFile Version 3.0" << endl;
  solution << "sol" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET STRUCTURED_POINTS" << endl;
  solution << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << endl;
  solution << "ORIGIN " << xmin << " " << ymin << " " << 0 << endl;
  solution << "SPACING " << hx << " " << hy << " " << 1 << endl;;
  solution << "POINT_DATA " << Nx*Ny << endl;
  solution << "SCALARS sol float" << endl;
  solution << "LOOKUP_TABLE default" << endl;
  for(int j=0; j<Ny; ++j)
    {
      for(int i=0; i<Nx; ++i)
	{
	  solution << _sol(i+j*Nx) << " ";
	}
      solution << endl;
    }
  solution.close();
}


VectorXd TimeScheme::GetSolution(){
  return _sol;
}


void EulerScheme::Integrate()
{
  _lab->BuildF(_sol,_f);
  _sol += _df->Get_dt()*_lab->Get_p();
  _t += _dt;
  _n+=1;
  }



#define _TIME_SCHEME_CPP
#endif
