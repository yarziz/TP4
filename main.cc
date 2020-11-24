#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Please, enter the name of your data file." << endl;
      exit(0);
    }
  // ----------------------- Fichier de données --------------------------------
  DataFile* df = new DataFile(argv[1]);
  df->ReadDataFile();

  // ---------------------------- Résolution  ----------------------------------
  // Pointeur contenant toutes les fonctions utiles
  Function* function = new Function(df);
  // Pointeur vers la classe Laplacian (discrétisation en espace)
  Laplacian* lap = new Laplacian(function, df);
  // Pointeur vers la classe TimeScheme (discrétisation en temps=)
  TimeScheme* time_scheme;

  int int_euler_scheme(1);
  if (df->Get_scheme() == "ImplicitEuler") {int_euler_scheme = 2;}

  switch(int_euler_scheme)
    {
    case 1:
      time_scheme = new EulerScheme(df,lap);
      break;
    case 2:
      //time_scheme = new ImplicitEulerScheme(df,lap);
      break;
    default:
      cout << "Ce choix n'est pas possible ! Veuillez recommencer !" << endl;
      exit(0);
    }

  time_scheme->Initialise();
  time_scheme->SaveSolution(0);

  int nb_iterations = int(ceil((df->Get_tfinal()-df->Get_t0())/df->Get_dt()));
  for (int n = 1; n < nb_iterations+1; n++)
    { // Boucle en temps
      time_scheme->Integrate();
      time_scheme->SaveSolution(n);
    }
  
  VectorXd exact_sol = lap->ExactSolution(df->Get_tfinal());
  VectorXd approx_sol = time_scheme->GetSolution();
  double error = ((approx_sol-exact_sol).array().abs()).maxCoeff();
  cout << "Erreur = " << error << endl;




  // ----------------------- Libère la mémoire ---------------------------------
  // Ne pas oublier de détruire tous les pointeurs que nous avons utilisé
  delete df; delete function; delete lap; delete time_scheme;

  return 0;
}
