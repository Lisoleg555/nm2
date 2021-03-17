#include "lib.h"

int main() {
	double max = M_PI, min = 0, maxT = 0;
	cout << "Input max T: ";
	cin >> maxT;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	int N, K, M;
	cout << "Input max number of iterations for h_x: ";
	cin >> N;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	cout << "Input max number of iterations for h_y: ";
	cin >> K;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	cout << "Input max number of iterations for tau: ";
	cin >> M;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	
	double hx = (double) max / N, hy = (double) max / K, tau = (double) maxT / M;
	vector<vector<vector<double>>> ukj(M + 1), correct(M + 1);
	for (int i = 0; i < M + 1; ++i) {
		ukj[i].resize(K + 1);
		correct[i].resize(K + 1);
		for (int j = 0; j < K + 1; ++j) {
			ukj[i][j].resize(N + 1, 0);
			correct[i][j].resize(N + 1, 0);
		}
	}
	vector<vector<double>> xY(3);
	xY[0].resize(N + 1, 0);
	xY[1].resize(K + 1, 0);
	xY[2].resize(M + 1, 0);
	for (int i = 1; i < N + 1; ++i)
		xY[0][i] = min + hx * i;
	for (int i = 1; i < K + 1; ++i)
		xY[1][i] = min + hy * i;
	for (int i = 1; i < M + 1; ++i)
		xY[2][i] = min + tau * i;
	xY[0].back() = max;
	xY[1].back() = max;
	xY[2].back() = maxT;
	char input;
	double a = 1, b = 1, mu = 2;
	cout << "Select parameters:\n1)a = 1, b = 1, mu = 1.\n2)a = 2, b = 1, mu = 1.\n3)a = 1, b = 2, mu = 1.\nelse)a = 1, b = 1, mu = 2.\n";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == '1')
		mu = 1;
	else if (input == '2') {
		a = 2;
		mu = 1;
	}
	else if (input == '3') {
		b = 2;
		mu = 1;
	}

	cout << "Select count of time plots (start and end already included)\n1)Every 2nd\n2)Every 10th\n3)Center one\n4)All\n";
	cin >> K;
	if (!cin || K < 1 || K > 4) {
		cout << "Wrong data\n";
		return -1;
	}
	vector<int> time;
	CountTimePlots(time, xY[2].size(), K);
	
	cout << "Analitical solve:\n";
	Analitical(xY, correct, mu);
	Gnuplot(xY, correct, time, "Analitical solve for time = ");
	
	cout << "Do you want to see the table with function values for analitical solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xY, correct);

	cout << "Alternating direction solve:\n";
	AltDir(xY, ukj, a, b, mu, hx, hy, tau);
	Gnuplot(xY, ukj, time, "Alternating direction solve for time = ");
	
	cout << "Do you want to see the table with function values for\nalternating direction solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xY, ukj);

	cout << "Do you want to see the error of calculations? (y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		ErrorsCalc(xY, ukj, correct);
	cout << "Fractional step solve:\n";
	FracStep(xY, ukj, a, b, mu, hx, hy, tau);
	Gnuplot(xY, ukj, time, "Fractional step solve for time = ");
	cout << "Do you want to see the table with function values for\nfractional step solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xY, ukj);
		
	cout << "Do you want to see the error of calculations? (y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		ErrorsCalc(xY, ukj, correct);
	return 0;
}
