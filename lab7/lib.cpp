#include "lib.h"

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
	file.close();
	FILE* gp = popen("gnuplot -persist", "w");
	if (gp == NULL) {
		cout << "Error: gnuplot didn't open\n";
		exit(-1);
	}
	string xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = std::to_string(xT[0][0]);
	xmax = std::to_string(xT[0].back() + 0.1);
	zmin = std::to_string(min);
	zmax = std::to_string(max);
	string gnuplot = "set title \"" + name + "\"\nset key off\n" + "set xrange [" + xmin + ":" + xmax +"]\nset yrange [" + xmin + ":" + xmax + "]\nset zrange [" + zmin + ":" + zmax + "]\n" + "set ylabel \"y\"\nset zlabel \"u\"\nset xlabel \"x\"\nsplot 'plot.txt' using 1:2:3 pt 7 ps 1 palette\n";
	fprintf(gp, gnuplot.c_str());
	pclose(gp);
	return;
}

void ErrorsCalc(vector<vector<double>>& xT, vector<vector<double>>& ukj, vector<vector<double>>& correct) {
	double min = abs(correct[1][0] - ukj[1][0]), max = abs(correct[1][0] - ukj[1][0]);
	cout << "Absolute errors:\n";
	for (int i = 0; i < xT[1].size(); ++i)
		for (int j = 0; j < xT[0].size(); ++j) {
			cout << "x = " << xT[0][j] << " y = " << xT[1][i] << " abs error = "
				 << abs(correct[i][j] - ukj[i][j]) << "\n";
			if (abs(correct[i][j] - ukj[i][j]) < min)
				min = abs(correct[i][j] - ukj[i][j]);
			if (abs(correct[i][j] - ukj[i][j]) > max)
				max = abs(correct[i][j] - ukj[i][j]);
		}
	cout << "min abs error = " << min << "\nmax abs error = " << max << "\n";
	return;
}

void PrintTables(vector<vector<double>>& xT, vector<vector<double>>& ukj) {
	for (int i = 0; i < xT[1].size(); ++i) {
		cout << "y = " << xT[1][i] << ":\n";
		for (int j = 0; j < ukj[i].size(); ++j)
			cout << "\tx = " << xT[0][j] << "\tu = " << ukj[i][j] << "\n";
	}
	cout << "----------------------------------------------------------------------\n";
	return;
}

void Analitical(vector<vector<double>>& xT, vector<vector<double>>& ukj) {
	for (int i = 0; i < xT[1].size(); ++i)
		for (int j = 0; j < xT[0].size(); ++j)
			ukj[i][j] = exp(-xT[1][i] - xT[0][j]) * cos(xT[1][i]) * cos(xT[0][j]);
	return;
}

void PrepareU(vector<vector<double>>& xT, vector<vector<double>>& u) {
	u[0][0] = exp(-xT[0][0]) * cos(xT[0][0]);
	u.back()[0] = 0;
	u[0].back() = 0;
	u.back().back() = 0;
	for (int i = 1; i < xT[1].size() - 1; ++i) {
		u[i][0] = exp(-xT[1][i]) * cos(xT[1][i]);
		u[i].back() = 0;
	}
	for (int i = 1; i < xT[0].size() - 1; ++i) {
		u[0][i] = exp(-xT[0][i]) * cos(xT[0][i]);
		u.back()[i] = 0;
	}
	for (int i = 1; i < xT[1].size() - 1; ++i) 
		for (int j = 1; j < xT[0].size() - 1; ++j)
			u[i][j] = u[i][0] * (xT[0].back() - xT[0][j]) / xT[0].back();
	return;
}

double FindMax(vector<vector<double>>& u) {
	double max = abs(u[1][1]);
	for (int i = 1; i < u.size(); ++i)
		for (int j = 1; j < u[i].size(); ++j)
			if (max < abs(u[i][j]))
				max = abs(u[i][j]);
	return max;
}

void SimpleIterations(vector<vector<double>>& u, double h, double tau, double eps, int kMax) {
	double max = 0, tau2 = pow(tau, 2), h2 = pow(h, 2);
	int k = 0;
	vector<vector<double>> uk = u;
	do {
		max = FindMax(u);
		for (int i = 1; i < u.size() - 1; ++i)
			for (int j = 1; j < u[i].size() - 1; ++j)
				uk[i][j] = (tau2 * u[i][j - 1] * (h - 1) - tau2 * u[i][j + 1] * (h + 1) - h2 * u[i + 1][j] *
						   (tau + 1) + h2 * u[i - 1][j] * (tau - 1)) / (4 * h2 * tau2 - 2 * tau2 - 2 * h2);
		u = uk;
		++k;
	} while(abs(max - FindMax(uk)) > eps && k < kMax);
	cout << k << " simple iteratons were made\n";
	return;
}

void Seidel(vector<vector<double>>& u, double h, double tau, double eps, int kMax) {
	double max = 0, tau2 = pow(tau, 2), h2 = pow(h, 2);
	int k = 0;
	vector<vector<double>> uk = u;
	do {
		max = FindMax(u);
		for (int i = 1; i < u.size() - 1; ++i)
			for (int j = 1; j < u[i].size() - 1; ++j)
				uk[i][j] = (tau2 * uk[i][j - 1] * (h - 1) - tau2 * u[i][j + 1] * (h + 1) - h2 * u[i + 1][j] *
						   (tau + 1) + h2 * uk[i - 1][j] * (tau - 1)) / (4 * h2 * tau2 - 2 * tau2 - 2 * h2);
		u = uk;
		++k;
	} while(abs(max - FindMax(uk)) > eps && k < kMax);
	cout << k << " seidel iteratons were made\n";
	return;
}

void SimpleIterationsHigh(vector<vector<double>>& u, double h, double tau, double eps, int kMax, double w) {
	double max = 0, tau2 = pow(tau, 2), h2 = pow(h, 2);
	int k = 0;
	vector<vector<double>> uk = u;
	do {
		max = FindMax(u);
		for (int i = 1; i < u.size() - 1; ++i)
			for (int j = 1; j < u[i].size() - 1; ++j) {
				uk[i][j] = (tau2 * uk[i][j - 1] * (h - 1) - tau2 * u[i][j + 1] * (h + 1) - h2 * u[i + 1][j] *
						   (tau + 1) + h2 * uk[i - 1][j] * (tau - 1)) / (4 * h2 * tau2 - 2 * tau2 - 2 * h2);
				uk[i][j] = uk[i][j] * w - (1 - w) * u[i][j];
			}
		u = uk;
		++k;
	} while(abs(max - FindMax(uk)) > eps && k < kMax);
	cout << k << " simple iteratons with high relaxation were made\n";
	return;
}
