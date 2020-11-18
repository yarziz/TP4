#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>

class DataFile {
private:
  // Tous les paramètres contenus dans le fichier de données
  double _xmin;
  double _xmax;
  double _hx;
  int _Nx;
  double _ymin;
  double _ymax;
  double _hy;
  int _Ny;
  double _sigma;
  double _t0;
  double _tfinal;
  double _dt;
  std::string _scheme;
  std::string _results;
  std::string _LBC, _RBC, _DBC, _UBC;
  // Nom du fichier de données
  const std::string _file_name;
  // Pour savoir si l'utilisateur a donné la valeur du paramètres
  // ou les paramètres par défaut doivent être utilisés
  bool _if_xmin;
  bool _if_xmax;
  bool _if_hx;
  bool _if_ymin;
  bool _if_ymax;
  bool _if_hy;
  bool _if_sigma;
  bool _if_t0;
  bool _if_tfinal;
  bool _if_dt;
  bool _if_scheme;
  bool _if_results;
  bool _if_LBC, _if_RBC, _if_DBC, _if_UBC;

public:
  // Constructeur
  DataFile(std::string file_name);
  // La lecture du fichier
  void ReadDataFile();
  // Nettoyer une ligne du fichier
  std::string clean_line(std::string &s);
  // Toutes les fonctions pour renvoyer les paramètres
  // Comme ce sont juste des fonctions Get, il est plus simple de les
  // définir directement dans le .cpp
  const double & Get_xmin() const {return _xmin;};
  const double & Get_xmax() const {return _xmax;};
  const double & Get_hx() const {return _hx;};
  const int & Get_Nx() const {return _Nx;};
  const double & Get_ymin() const {return _ymin;};
  const double & Get_ymax() const {return _ymax;};
  const double & Get_hy() const {return _hy;};
  const int & Get_Ny() const {return _Ny;};
  const double & Get_sigma() const {return _sigma;};
  const double & Get_t0() const {return _t0;};
  const double & Get_tfinal() const {return _tfinal;};
  const double & Get_dt() const {return _dt;};
  const std::string & Get_scheme() const {return _scheme;};
  const std::string & Get_results() const {return _results;};
  const std::string & Get_LBC() const {return _LBC;};
  const std::string & Get_RBC() const {return _RBC;};
  const std::string & Get_DBC() const {return _DBC;};
  const std::string & Get_UBC() const {return _UBC;};
};

#define _DATA_FILE_H
#endif
