#include "../common/scatter_matrix.h"

double Phi(double, double, double, double);
void Gnuplot(vector<vector<double>>&, vector<vector<double>>&, string);
void ErrorsCalc(vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&);
void PrintTables(vector<vector<double>>&, vector<vector<double>>&);
double InitialConditions(double);
void Analitical(vector<vector<double>>&, vector<vector<double>>&, double, double, double);
void Explit(vector<vector<double>>&, vector<vector<double>>&, double, double, double, double, double, bool);
void Implict(vector<vector<double>>&, vector<vector<double>>&, double, double, double, double, double, bool);
void CrankNicholson(vector<vector<double>>&, vector<vector<double>>&, double, double, double, double, double, bool);
