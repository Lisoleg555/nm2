#include "../common/scatter_matrix.h"

double F(double, double, double, double, double, double);
void CountTimePlots(vector<int>&, int, int);
void Gnuplot(vector<vector<double>>&, vector<vector<vector<double>>>&, vector<int>&, string);
void ErrorsCalc(vector<vector<double>>&, vector<vector<vector<double>>>&, vector<vector<vector<double>>>&);
void PrintTables(vector<vector<double>>&, vector<vector<vector<double>>>&);
void Analitical(vector<vector<double>>&, vector<vector<vector<double>>>&, double);
void AltDir(vector<vector<double>>&, vector<vector<vector<double>>>&, double, double, double, double, double, double);
void FracStep(vector<vector<double>>&, vector<vector<vector<double>>>&, double, double, double, double, double, double);
