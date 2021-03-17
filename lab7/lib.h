#include "../common/scatter_matrix.h"

void Gnuplot(vector<vector<double>>&, vector<vector<double>>&, string);
void ErrorsCalc(vector<vector<double>>&, vector<vector<double>>&, vector<vector<double>>&);
void PrintTables(vector<vector<double>>&, vector<vector<double>>&);
void PrepareU(vector<vector<double>>&, vector<vector<double>>&);
void Analitical(vector<vector<double>>&, vector<vector<double>>&);
void SimpleIterations(vector<vector<double>>&, double, double, double, int);
void Seidel(vector<vector<double>>&, double, double, double, int);
void SimpleIterationsHigh(vector<vector<double>>&, double, double, double, int, double);
