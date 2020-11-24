#ifndef _TIME_SCHEME_H

#include "Laplacian.h"

class TimeScheme
{

public:
  // Constructeur par défaut
  TimeScheme(DataFile* data_file, Laplacian* lap);
  // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
  // alors le compilateur en génère un implicitement.
  virtual ~TimeScheme();
  void SaveSolution(int n);
  void Initialise();
  Eigen::VectorXd GetSolution();
  virtual void Integrate()=0;

 protected:
  // Pointeur vers le fichier de données
  DataFile* _df;
  // Pointeur vers le laplacien
  Laplacian* _lap;

  int _n;
  double _t;
  Eigen::VectorXd _sol;
  Eigen::VectorXd _f;
};

class EulerScheme : public TimeScheme
{
public:
  EulerScheme(DataFile* data_file, Laplacian* lap);
  void Integrate();
};

class ImplicitEulerScheme : public TimeScheme
{
public:
  ImplicitEulerScheme(DataFile* data_file, Laplacian* lap);
};

#define _TIME_SCHEME_H
#endif
