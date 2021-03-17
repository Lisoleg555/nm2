#include "lib.h"

double Phi(double t, double a, double b, double c) {
	return exp((c - a) * t) * (cos(b * t) + sin(b * t));
}

void Gnuplot(vector<vector<double>>& xT, vector<vector<double>>& ukj, string name) {
	double min = ukj[0][0], max = ukj[0][0];
	for (int i = 0; i < ukj.size(); ++i)
		for (int j = 0; j < ukj[i].size(); ++j) {
			if (ukj[i][j] > max)
				max = ukj[i][j];
			if (ukj[i][j] < min)
				min = ukj[i][j];
		}
	ofstream file("plot.txt");
	for (int i = 0; i < xT[1].size(); ++i)
		for (int j = 0; j < xT[0].size(); ++j)
			file << xT[0][j] << " " << xT[1][i] << " " << ukj[i][j] << "\n";
	file << 0 << " " << 0 << " " << -0.01 << "\n";
	file.close();
	FILE* gp = popen("gnuplot -persist", "w");
	if (gp == NULL) {
		cout << "Error: gnuplot didn't open\n";
		exit(-1);
	}
	string xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = std::to_string(xT[0][0]);
	xmax = std::to_string(xT[0].back());
	ymin = std::to_string(xT[1][0]);
	ymax = std::to_string(xT[1].back());
	zmin = std::to_string(min - 0.1);
	zmax = std::to_string(max);
	string gnuplot = "set title \"" + name + "\"\nset key off\n" + "set xrange [" + xmin + ":" + xmax +"]\nset yrange [" + ymin + ":" + ymax + "]\nset zrange [" + zminn + ":" + zmax + "]\n" + "set ylabel \"t\"\nset zlabel \"u\"\nset xlabel \"x\"\nsplot 'plot.txt' using 1:2:3 pt 7 ps 1 palette\n";
	fprintf(gp, gnuplot.c_str());
	pclose(gp);
	return;
}

void ErrorsCalc(vector<vector<double>>& xT, vector<vector<double>>& ukj, vector<vector<double>>& correct) {
	double min = abs(abs(correct[1][0]) - abs(ukj[1][0])), max = -1;
	cout << "Absolute errors:\n";
	for (int i = 0; i < xT[1].size(); ++i)
		for (int j = 0; j < xT[0].size(); ++j) {
			cout << "x = " << xT[0][j] << " t = " << xT[1][i] << " abs error = "
				 << abs(abs(correct[i][j]) - abs(ukj[i][j])) << "\n";
			if (abs(abs(correct[i][j]) - abs(ukj[i][j])) < min)
				min = abs(abs(correct[i][j]) - abs(ukj[i][j]));
			if (abs(abs(correct[i][j]) - abs(ukj[i][j])) > max)
				max = abs(abs(correct[i][j]) - abs(ukj[i][j]));
		}
	cout << "min abs error = " << min << "\nmax abs error = " << max << "\n";
	return;
}

void PrintTables(vector<vector<double>>& xT, vector<vector<double>>& ukj) {
	for (int i = 0; i < xT[1].size(); ++i) {
		cout << "t = " << xT[1][i] << ":\n";
		for (int j = 0; j < ukj[i].size(); ++j)
			cout << "\tx = " << xT[0][j] << "\tu = " << ukj[i][j] << "\n";
	}
	cout << "----------------------------------------------------------------------\n";
	return;
}

double InitialConditions(double x) {
	return sin(x);
}

void Analitical(vector<vector<double>>& xT, vector<vector<double>>& ukj, double a, double b, double c) {
	for (int i = 0; i < xT[1].size(); ++i) 
		for (int j = 0; j < xT[0].size(); ++j)
			ukj[i][j] = exp((c - a) * xT[1][i]) * sin(xT[0][j] + b * xT[1][i]);
	return;
}

void Explit(vector<vector<double>>& xT, vector<vector<double>>& u, double a, double b, double c,
			double h, double tau, bool approx) {
	for (int i = 0; i < xT[1].size() - 1; ++i) {
		for (int j = 1; j < xT[0].size() - 1; ++j) {
			u[i + 1][j] = tau * (a * (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / pow(h, 2) +
						  b * (u[i][j + 1] - u[i][j - 1]) / 2 / h + c * u[i][j]) + u[i][j];
			
		}
		if (approx) {
			u[i + 1][0] = (h * Phi(xT[1][i], a, b, c) - u[i + 1][1]) / (h - 1);
			u[i + 1].back() = (u[i + 1][xT[0].size() - 2] - h * Phi(xT[1][i], a, b, c)) / (h + 1);
		}
		else {
			u[i + 1][0] = (2 * h * Phi(xT[1][i], a, b, c) + u[i + 1][2] - 4 * u[i + 1][1]) / (2 * h - 3);
			u[i + 1].back() = -(2 * h * Phi(xT[1][i], a, b, c) + u[i + 1][xT[0].size() - 3] -
							  4 * u[i + 1][xT[0].size() - 2]) / (2 * h + 3);
		}
	}
	return;
}

void Implict(vector<vector<double>>& xT, vector<vector<double>>& ukj, double a, double b, double c,
			 double h, double tau, bool approx) {
	vector<vector<double>> mat(xT[0].size());
	double del1 = 0, del2 = 0;
	for (int i = 1; i < mat.size() - 1; ++i)
		mat[i].resize(3);
	mat[0].resize(2);
	mat.back().resize(2);
	for (int i = 1; i < xT[0].size() - 1; ++i) {
		mat[i][2] = (a / pow(h, 2) + b / 2 / h) * tau;
		mat[i][1] = c * tau - 2 * a * tau / pow(h, 2) - 1;
		mat[i][0] = (a / pow(h, 2) - b * tau / 2 / h) * tau;
	}
	if (approx) {
		mat[0][0] = h - 1;
		mat[0][1] = 1;
		mat.back()[0] = -1;
		mat.back()[1] = h + 1;
	}
	else {
		del1 = -1 / mat[1][2];
		del2 = 1 / mat[mat.size() - 2][0];
		mat[0][0] = 2 * h - 3 - del1 * mat[1][0];
		mat[0][1] = 4 - del1 * mat[1][1];
		mat.back()[0] = -4 - del2 * mat[mat.size() - 2][1];
		mat.back()[1] = 2 * h + 3 - del2 * mat[mat.size() - 2][2];
	}
	
	vector<double> vec(xT[0].size());
	for (int i = 1; i < xT[1].size(); ++i) {
		for (int j = 1; j < vec.size() - 1; ++j)
			vec[j] = -ukj[i - 1][j];
		if (approx) {
			vec[0] = h * Phi(xT[1][i], a, b, c);
			vec.back() = -h * Phi(xT[1][i], a, b, c);
		}
		else {
			vec[0] = 2 * h * Phi(xT[1][i], a, b, c) - del1 * vec[1];
			vec.back() = -2 * h * Phi(xT[1][i], a, b, c) - del2 * vec[vec.size() - 2];
		}
		
		ThreeDiag help;
		help.PutMatrix(mat, vec);
		help.Progonka();
		ukj[i] = help.GetAnswer();
	}
	return;
}

void CrankNicholson(vector<vector<double>>& xT, vector<vector<double>>& u, double a, double b, double c,
					double h, double tau, bool approx) {
	vector<vector<double>> mat(xT[0].size());
	double del1 = 0, del2 = 0;
	for (int i = 1; i < mat.size() - 1; ++i)
		mat[i].resize(3);
	mat[0].resize(2);
	mat.back().resize(2);
	double tmp;
	for (int i = 1; i < xT[0].size() - 1; ++i) {
		mat[i][2] = -(a / pow(h, 2) + b / 2 / h) / 2;
		mat[i][1] = 1 / tau + a / pow(h, 2) + c / 2;
		mat[i][0] = (- a / pow(h, 2) + b / 2 / h) / 2;
	}
	if (approx) {
		mat[0][0] = h - 1;
		mat[0][1] = 1;
		mat.back()[0] = -1;
		mat.back()[1] = h + 1;
	}
	else {
		del1 = -1 / mat[1][2];
		del2 = 1 / mat[mat.size() - 2][0];
		mat[0][0] = 2 * h - 3 - del1 * mat[1][0];
		mat[0][1] = 4 - del1 * mat[1][1];
		mat.back()[0] = -4 - del2 * mat[mat.size() - 2][1];
		mat.back()[1] = 2 * h + 3 - del2 * mat[mat.size() - 2][2];
	}
	
	vector<double> vec(xT[0].size());
	for (int i = 1; i < xT[1].size(); ++i) {
		for (int j = 1; j < vec.size() - 1; ++j) {
			vec[j] = u[i - 1][j - 1] * (a / pow(h , 2) - b / h / 2) / 2 + u[i - 1][j] * (1 / tau -
					 a / pow(h, 2) + c / 2) + u[i - 1][j + 1] * (a / pow(h, 2) + b / 2 / h) / 2;
		}
		if (approx) {
			vec[0] = h * Phi(xT[1][i], a, b, c);
			vec.back() = -h * Phi(xT[1][i], a, b, c);
		}
		else {
			vec[0] = 2 * h * Phi(xT[1][i], a, b, c) - del1 * vec[1];
			vec.back() = -2 * h * Phi(xT[1][i], a, b, c) - del2 * vec[vec.size() - 2];
		}
		
		ThreeDiag help;
		help.PutMatrix(mat, vec);
		help.Progonka();
		u[i] = help.GetAnswer();
	}
	return;
}
