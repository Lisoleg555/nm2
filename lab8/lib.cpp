#include "lib.h"

double F(double x, double y, double t, double a, double b, double mu) {
	return sin(x) * sin(y) * (mu * cos(mu * t) + (a + b) * sin(mu * t));
}

void CountTimePlots(vector<int>& time, int tSize, int k) {
	switch(k) {
		case 1: {
			time.resize(tSize / 2 + 1);
			for (int i = 0; i < time.size() - 1; ++i)
				time[i] = i * 2;
			time.back() = tSize - 1;
			break;
		}
		case 2: {
			k = (tSize - 1) % 10;
			time.resize(tSize / 10 + (k == 0 ? 1 : 2));
			for (int i = 0; i < time.size() - 1; ++i)
				time[i] = i * 10;
			time.back() = tSize - 1;
			break;
		}
		case 3: {
			time.resize(3);
			time[0] = 0;
			time[1] = tSize / 2;
			time[2] = tSize - 1;
			break;
		}
		case 4: {
			time.resize(tSize);
			for (int i = 0; i < tSize; ++i)
				time[i] = i;
			break;
		}
	}
	return;
}

void Gnuplot(vector<vector<double>>& xT, vector<vector<vector<double>>>& ukj, vector<int>& time, string name) {
	for (int k = 0; k < time.size(); ++k) {
		double min = ukj[time[k]][0][0], max = ukj[k][0][0];
		for (int i = 0; i < ukj[time[k]].size(); ++i)
			for (int j = 0; j < ukj[time[k]][i].size(); ++j) {
				if (ukj[time[k]][i][j] > max)
					max = ukj[time[k]][i][j];
				if (ukj[k][i][j] < min)
					min = ukj[time[k]][i][j];
			}
		ofstream file("plot.txt");
		for (int i = 0; i < xT[1].size(); ++i)
			for (int j = 0; j < xT[0].size(); ++j)
				file << xT[0][j] << " " << xT[1][i] << " " << ukj[time[k]][i][j] << "\n";
		file.close();
		FILE* gp = popen("gnuplot -persist", "w");
		if (gp == NULL) {
			cout << "Error: gnuplot didn't open\n";
			exit(-1);
		}
		string xmin, xmax, ymin, ymax, zmin, zmax, timeName;
		timeName = std::to_string(xT[2][time[k]]);
		xmin = std::to_string(xT[0][0]);
		xmax = std::to_string(xT[0].back());
		zmin = std::to_string(min - 0.1);
		zmax = std::to_string(max + 0.1);
		string gnuplot = "set title \"" + name + timeName + "\"\nset key off\n" + "set xrange [" + xmin + ":" + xmax +"]\nset yrange [" + xmin + ":" + xmax + "]\nset zrange [" + zmin + ":" + zmax + "]\n" + "set ylabel \"y\"\nset zlabel \"u\"\nset xlabel \"x\"\nsplot 'plot.txt' using 1:2:3 pt 7 ps 1 palette\n";
		fprintf(gp, gnuplot.c_str());
		pclose(gp);
	}
	return;
}

void ErrorsCalc(vector<vector<double>>& xT, vector<vector<vector<double>>>& ukj,
				vector<vector<vector<double>>>& correct) {
	double min = abs(correct[1][1][0] - ukj[1][1][0]), max = abs(correct[1][1][0] - ukj[1][1][0]);
	cout << "Absolute errors:\n";
	for (int k = 0; k < xT[2].size(); ++k)
		for (int i = 0; i < xT[1].size(); ++i)
			for (int j = 0; j < xT[0].size(); ++j) {
				cout << "x = " << xT[0][j] << " y = " << xT[1][i] << " t = " << xT[2][k] << " abs error = "
					 << abs(correct[k][i][j] - ukj[k][i][j]) << "\n";
				if (abs(correct[k][i][j] - ukj[k][i][j]) < min)
					min = abs(correct[k][i][j] - ukj[k][i][j]);
				if (abs(correct[k][i][j] - ukj[k][i][j]) > max)
					max = abs(correct[k][i][j] - ukj[k][i][j]);
			}
	cout << "min abs error = " << min << "\nmax abs error = " << max << "\n";
	return;
}

void PrintTables(vector<vector<double>>& xT, vector<vector<vector<double>>>& ukj) {
	for (int k = 0; k < xT[2].size(); ++k) {
		cout << "t = " << xT[2][k] << ":\n";
		for (int i = 0; i < xT[1].size(); ++i) {
			cout << "\ty = " << xT[1][i] << ":\n";
			for (int j = 0; j < xT[0].size(); ++j)
				cout << "\t\tx = " << xT[0][j] << "\tu = " << ukj[k][i][j] << "\n";
		}
	}
	cout << "----------------------------------------------------------------------\n";
	return;
}

void Analitical(vector<vector<double>>& xT, vector<vector<vector<double>>>& ukj, double mu) {
	for (int k = 0; k < xT[2].size(); ++k)
		for (int i = 0; i < xT[1].size(); ++i)
			for (int j = 0; j < xT[0].size(); ++j)
				ukj[k][i][j] = sin(xT[0][j]) * sin(xT[1][i]) * sin(mu * xT[2][k]);
	return;
}

void AltDir(vector<vector<double>>& xT, vector<vector<vector<double>>>& u, double a, double b, double mu,
			double hx, double hy, double tau) {
	double tmpY = b * tau / pow(hy, 2) / 2;
	double tmpX = a * tau / pow(hx, 2) / 2;
	double tau2 = tau / 2;
	vector<vector<double>> matx(xT[0].size() - 2);
	for (int i = 1; i < matx.size() - 1; ++i)
		matx[i].resize(3);
	matx[0].resize(2);
	matx.back().resize(2);
	matx[1][2] = -tmpX;
	matx[1][1] = 2 * tmpX + 1;
	matx[1][0] = matx[1][2];
	for (int i = 2; i < matx.size() - 1; ++i) {
		matx[i][2] = matx[1][2];
		matx[i][1] = matx[1][1];
		matx[i][0] = matx[1][0];
	}
	matx[0][0] = matx[1][1];
	matx[0][1] = matx[1][2];
	matx.back()[0] = matx[1][0];
	matx.back()[1] = matx[1][1];

	vector<vector<double>> maty(xT[1].size() - 2);
	for (int i = 1; i < maty.size() - 1; ++i)
		maty[i].resize(3);
	maty[0].resize(2);
	maty.back().resize(2);
	maty[1][2] = -tmpY;
	maty[1][1] = 2 * tmpY + 1;
	maty[1][0] = maty[1][2];
	for (int i = 2; i < maty.size() - 1; ++i) {
		maty[i][2] = maty[1][2];
		maty[i][1] = maty[1][1];
		maty[i][0] = maty[1][0];
	}
	maty[0][0] = maty[1][1];
	maty[0][1] = maty[1][2];
	maty.back()[0] = maty[1][0];
	maty.back()[1] = maty[1][1];
	vector<vector<double>> ukj(xT[1].size() - 2);
	vector<double> vecx(xT[0].size() - 2);
	vector<double> vecy(xT[1].size() - 2);
	for (int k = 0; k < xT[2].size() - 1; ++k) {
		//прогонка для к+1/2
		for (int i = 1; i <= ukj.size(); ++i) {
			for (int j = 1; j <= vecx.size(); ++j)
				vecx[j - 1] = u[k][i][j] + tmpY * (u[k][i + 1][j] - 2 * u[k][i][j] + u[k][i - 1][j]) +
							  tau2 * F(xT[0][j], xT[1][i], xT[2][k] + tau2, a, b, mu);
			ThreeDiag help;
			help.PutMatrix(matx, vecx);
			help.Progonka();
			ukj[i - 1] = help.GetAnswer();
		}
		//прогонка для к + 1
		
		for (int i = 0; i < vecy.size(); ++i)
			vecy[i] = ukj[i][0] +  tmpX * (ukj[i][1] - 2 * ukj[i][0]) +
					  tau2 * F(xT[0][1], xT[1][i], xT[2][k] + tau2, a, b, mu);
		ThreeDiag help1;
		help1.PutMatrix(maty, vecy);
		help1.Progonka();
		vecy = help1.GetAnswer();
		for (int i = 0; i < vecy.size(); ++i)
			u[k + 1][i + 1][1] = vecy[i];
		for (int i = 1; i < xT[0].size() - 2; ++i) {
			for (int j = 0; j < vecy.size(); ++j)
				vecy[j] = ukj[j][i] + tmpX * (ukj[j][i + 1] - 2 * ukj[j][i] + ukj[j][i - 1]) +
						  tau2 * F(xT[0][i + 1], xT[1][j + 1], xT[2][k + 1] + tau2, a, b, mu);
			ThreeDiag help2;
			help2.PutMatrix(maty, vecy);
			help2.Progonka();
			vecy = help2.GetAnswer();
			for (int j = 0; j < vecy.size(); ++j)
				u[k + 1][j + 1][i + 1] = vecy[j];
		}
	}
	return;
}

void FracStep(vector<vector<double>>& xT, vector<vector<vector<double>>>& u, double a, double b, double mu,
			double hx, double hy, double tau) {
	double tmpY = b * tau / pow(hy, 2);
	double tmpX = a * tau / pow(hx, 2);
	double tau2 = tau / 2;
	vector<vector<double>> matx(xT[0].size() - 2);
	for (int i = 1; i < matx.size() - 1; ++i)
		matx[i].resize(3);
	matx[0].resize(2);
	matx.back().resize(2);
	matx[1][2] = -tmpX;
	matx[1][1] = 2 * tmpX + 1;
	matx[1][0] = matx[1][2];
	for (int i = 2; i < matx.size() - 1; ++i) {
		matx[i][2] = matx[1][2];
		matx[i][1] = matx[1][1];
		matx[i][0] = matx[1][0];
	}
	matx[0][0] = matx[1][1];
	matx[0][1] = matx[1][2];
	matx.back()[0] = matx[1][0];
	matx.back()[1] = matx[1][1];

	vector<vector<double>> maty(xT[1].size() - 2);
	for (int i = 1; i < maty.size() - 1; ++i)
		maty[i].resize(3);
	maty[0].resize(2);
	maty.back().resize(2);
	maty[1][2] = -tmpY;
	maty[1][1] = 2 * tmpY + 1;
	maty[1][0] = maty[1][2];
	for (int i = 2; i < maty.size() - 1; ++i) {
		maty[i][2] = maty[1][2];
		maty[i][1] = maty[1][1];
		maty[i][0] = maty[1][0];
	}
	maty[0][0] = maty[1][1];
	maty[0][1] = maty[1][2];
	maty.back()[0] = maty[1][0];
	maty.back()[1] = maty[1][1];
	vector<vector<double>> ukj(xT[1].size() - 2);
	vector<double> vecx(xT[0].size() - 2);
	vector<double> vecy(xT[1].size() - 2);
	for (int k = 0; k < xT[2].size() - 1; ++k) {
		//прогонка для к+1/2
		for (int i = 1; i <= ukj.size(); ++i) {
			for (int j = 1; j <= vecx.size(); ++j)
				vecx[j - 1] = u[k][i][j] + tau2 * F(xT[0][j], xT[1][i], xT[2][k], a, b, mu);
			ThreeDiag help;
			help.PutMatrix(matx, vecx);
			help.Progonka();
			ukj[i - 1] = help.GetAnswer();
		}
		//прогонка для к + 1
		
		for (int i = 0; i < vecy.size(); ++i)
			vecy[i] = ukj[i][0] + tau2 * F(xT[0][1], xT[1][i + 1], xT[2][k + 1], a, b, mu);
		ThreeDiag help1;
		help1.PutMatrix(maty, vecy);
		help1.Progonka();
		vecy = help1.GetAnswer();
		for (int i = 0; i < vecy.size(); ++i)
			u[k + 1][i + 1][1] = vecy[i];
		for (int i = 1; i < xT[0].size() - 1; ++i) {
			for (int j = 0; j < vecy.size(); ++j)
				vecy[j] = ukj[j][i] + tau2 * F(xT[0][i + 1], xT[1][j + 1], xT[2][k + 1], a, b, mu);
			ThreeDiag help2;
			help2.PutMatrix(maty, vecy);
			help2.Progonka();
			vecy = help2.GetAnswer();
			for (int j = 0; j < vecy.size(); ++j)
				u[k + 1][j + 1][i + 1] = vecy[j];
		}
	}
	return;
}
