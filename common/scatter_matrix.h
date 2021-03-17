#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#define _USE_MATH_DEFINES

using namespace std;

const std::string zminn = std::to_string(-0.2);

class ThreeDiag {
public:
	ThreeDiag();
	void PutMatrix(vector<vector<double>>&, vector<double>&);
	void GetAll();
	void GetFromFile(string);
	void GetVector();
	void PrepareSpline(vector<double>&, vector<double>&);
	void GetMatrix();
	void Progonka();
	void CheckSolve(vector<double>&);
	void PutEndDiff(vector<vector<double>>&);
	vector<double>& GetAnswer();
	void PrintAnswer();
	void PrintSystem();
	void PrintSpline();
	double Det();
	vector<double> operator [] (int) const;
	friend ostream& operator << (ostream&, const ThreeDiag&);
	int KnowSize();
	~ThreeDiag();
private:
	int size;
	vector<double> answer;
	vector<vector<double>> elements;
};
