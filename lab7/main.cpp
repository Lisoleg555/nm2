#include "lib.h"

int main() {
	double max = M_PI / 2, min = 0;
	int N, K;
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
	double eps = 0;
	cout << "Input eps for iterations methods: ";
	cin >> eps;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	
	double h = (double) max / N, tau = (double) max / K;
	vector<vector<double>> ukj(K + 1), correct(K + 1);
	for (int i = 0; i < K + 1; ++i) {
		ukj[i].resize(N + 1, 0);
		correct[i].resize(N + 1, 0);
	}
	vector<vector<double>> xY(2);
	xY[0].resize(N + 1);
	xY[1].resize(K + 1);
	for (int i = 1; i < N + 1; ++i) {
		xY[0][i] = min + h * i;
		xY[0][i] = min + h * i;
	}
	for (int i = 1; i < K + 1; ++i) {
		xY[1][i] = min + tau * i;
		xY[1][i] = min + tau * i;
	}
	xY[0][0] = min;
	xY[1][0] = min;
	xY[0].back() = max;
	xY[1].back() = max;
	
	cout << "Analitical solve:\n";
	Analitical(xY, correct);
	Gnuplot(xY, correct, "Analitical solve");
	cout << "Do you want to see the table with function values for analitical solve(y\\n): ";
	char input;
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xY, correct);

	cout << "Input max number of iterations for algorithms: ";
	cin >> K;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}

	cout << "Simple iterations solve:\n";
	PrepareU(xY, ukj);
	SimpleIterations(ukj, h, tau, eps, K);
	Gnuplot(xY, ukj, "Simple iterations solve");

	cout << "Do you want to see the table with function values for\nsimple iterations solve(y\\n): ";
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

	cout << "Seidel solve:\n";
	PrepareU(xY, ukj);
	Seidel(ukj, h, tau, eps, K);
	Gnuplot(xY, ukj, "Seidel solve");
	cout << "Do you want to see the table with function values for Seidel solve(y\\n): ";
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

	cout << "Simple iterations with high relaxation solve:\n";
	PrepareU(xY, ukj);
	cout << "Input parameter for high relaxation: ";
	cin >> max;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	SimpleIterationsHigh(ukj, h, tau, eps, K, max);
	Gnuplot(xY, ukj, "Simple iterations with high relaxation solve");
	cout << "Do you want to see the table with function values for\nsimple iterations with high relaxation solve(y\\n): ";
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
