#include "../common/scatter_matrix.h"

double F(double, double);
double Phi(double);
void Gnuplot(vector<vector<double>>&, vector<vector<double>>&, string);
void ErrorsCalc(vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&);
void PrintTables(vector<vector<double>>&, vector<vector<double>>&);
double InitialConditions(double);
void Analitical(vector<vector<double>>&, vector<vector<double>>&);
void Explit(vector<vector<double>>&, vector<vector<double>>&, double, double, bool);
void Implict(vector<vector<double>>&, vector<vector<double>>&, double, double, bool);
