#include "lib.h"

int main() {
	double maxX = M_PI, minT = 0, minX = 0;
	double maxT;
	cout << "Input max T: ";
	cin >> maxT;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	int N, K;
	cout << "Input max number of iterations for h: ";
	cin >> N;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	cout << "Input max number of iterations for tau: ";
	cin >> K;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	double a, b, c;
	cout << "Input a > 0, b > 0 and c < 0:\n";
	cin >> a >> b >> c;
	if (a <= 0 || b <= 0 || c >= 0) {
		cout << "Wrong data\n";
		return -1;
	}
	bool approx = true;
	char input;
	cout << "Do you want to use approximation of second order input(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		approx = false;
	double h = (double) maxX / N, tau = (double) maxT / K;
	vector<vector<double>> ukj(K + 1), correct(K + 1);
	for (int i = 0; i < K + 1; ++i) {
		ukj[i].resize(N + 1, 0);
		correct[i].resize(N + 1, 0);
	}
	vector<vector<double>> xTau(2);
	xTau[0].resize(N + 1);
	xTau[1].resize(K + 1);
	for (int i = 1; i < N; ++i)
		xTau[0][i] = minX + h * i;
	for (int i = 1; i < K; ++i)
		xTau[1][i] = minT + tau * i;
	xTau[0][0] = minX;
	xTau[0].back() = maxX;
	xTau[1][0] = minT;
	xTau[1].back() = maxT;
	cout << "Analitical solve:\n";
	Analitical(xTau, correct, a, b, c);
	
	Gnuplot(xTau, correct, "Analitical solve");
	cout << "Do you want to see the table with function values for analitical solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xTau, correct);

	for (int i = 0; i < N + 1; ++i)
		ukj[0][i] = InitialConditions(xTau[0][i]);
	cout << "Explit solve:\n";
	Explit(xTau, ukj, a, b, c, h, tau, approx);
	Gnuplot(xTau, ukj, "Explit solve");
	cout << "Do you want to see the table with function values for explict solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xTau, ukj);

	cout << "Do you want to see the error of calculations? (y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		ErrorsCalc(xTau, ukj, correct);

	cout << "Implict solve:\n";
	Implict(xTau, ukj, a, b, c, h, tau, approx);
	Gnuplot(xTau, ukj, "Implict solve");
	cout << "Do you want to see the table with function values for implict solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xTau, ukj);
		
	cout << "Do you want to see the error of calculations? (y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		ErrorsCalc(xTau, ukj, correct);
	
	cout << "Crank Nicholson solve:\n";
	CrankNicholson(xTau, ukj, a, b, c, h, tau, approx);
	Gnuplot(xTau, ukj, "Crank Nicholson solve");
	cout << "Do you want to see the table with function values for Crank Nicholson solve(y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		PrintTables(xTau, ukj);
	cout << "Do you want to see the error of calculations? (y\\n): ";
	cin >> input;
	if (!cin) {
		cout << "Wrong data\n";
		return -1;
	}
	if (input == 'y')
		ErrorsCalc(xTau, ukj, correct);
	return 0;
}
