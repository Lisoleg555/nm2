#include "scatter_matrix.h"

ThreeDiag::ThreeDiag() {
	size = 0;
}

void ThreeDiag::PutMatrix(vector<vector<double>>& mat, vector<double>& ans) {
	this->size = mat.size();
	this->elements.resize(this->size);
	this->answer.resize(this->size);
	for (int i = 1; i < this->size - 1; ++i) {
		this->elements[i].resize(3);
		this->elements[i][0] = mat[i][0];
		this->elements[i][1] = mat[i][1];
		this->elements[i][2] = mat[i][2];
		this->answer[i] = ans[i];
	}
	this->elements[0].resize(2);
	this->elements.back().resize(2);
	this->elements[0][0] = mat[0][0];
	this->elements[0][1] = mat[0][1];
	this->elements.back()[0] = mat.back()[0];
	this->elements.back()[1] = mat.back()[1];
	this->answer[0] = ans[0];
	this->answer.back() = ans.back();
	return;
}

double ThreeDiag::Det() {
	double f = elements[0][0];
	double ff = 1;
	for (int i = 1; i < size; ++i) {
		double tmp = f * elements[i][1] - elements[i][0] * elements[i - 1].back() * ff;
		ff = f;
		f = tmp;
	}
	return f;
}

void ThreeDiag::GetAll() {
	cout << "input size\n";
	cin >> size;
	elements.resize(size);
	for (int i = 1; i < size - 1; ++i)
		elements[i].resize(3);
	elements[0].resize(2);
	elements.back().resize(2);
	answer.resize(size);
	cout << "input matrix and vector\n";
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < elements[i].size(); ++j)
			cin >> elements[i][j];
		cin >> answer[i];
	}
	return;
}

void ThreeDiag::GetFromFile(string fileName) {
	ifstream file;
	file.open(fileName);
	if (!file.is_open()) {
		cout << "file isn't open\n";
		return;
	}
	file >> size;
	elements.resize(size);
	for (int i = 1; i < size - 1; ++i)
		elements[i].resize(3);
	elements[0].resize(2);
	elements.back().resize(2);
	answer.resize(size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < elements[i].size(); ++j)
			file >> elements[i][j];
		file >> answer[i];
	}
	file.close();
	return;
}

void ThreeDiag::GetVector() {
	cout << "input vector with " << size << " elements\n";
	for (int i = 0; i < size; ++i)
		cin >> answer[i];
	return;
}

void ThreeDiag::PrepareSpline(vector<double>& h, vector<double>& f) {
	size = f.size() - 2;
	elements.resize(size);
	for (int i = 1; i < size - 1; ++i)
		elements[i].resize(3);
	elements[0].resize(2);
	elements.back().resize(2);
	answer.resize(size);
	elements[0][0] = 2 * (h[0] + h[1]);
	elements[0][1] = h[1];
	answer[0] = 3 * ((f[2] - f[1]) / h[1] - (f[1] - f[0]) / h[0]);
	for (int i = 1; i < size - 1; ++i) {
		elements[i][0] = h[i];
		elements[i][1] = 2 * (h[i + 1] + h[i]);
		elements[i][2] = h[i + 1];
		answer[i] = 3 * ((f[i + 2] - f[i + 1]) / h[i + 1] - (f[i + 1] - f[i]) / h[i]);
	}
	elements.back()[0] = h[size - 1];
	elements.back()[1] = 2 * (h[size] + h[size - 1]);
	answer.back() = 3 * ((f.back() - f[size]) / h[size] - (f[size] - f[size - 1]) / h[size - 1]);
	return;
}

void ThreeDiag::PrintSpline() {
	int number = 1;
	if (elements[0][0] != 0) {
		if (elements[0][0] == 1)
			cout << "c2 ";
		else if (elements[0][0] == -1)
			cout << "-c2 ";
		else
			cout << elements[0][0] << "c2 ";
		if (elements[0][1] == 1)
			cout << "+ c3 ";
		else if (elements[0][1] == -1)
			cout << "- c3 ";
		else if (elements[0][1] < 0)
			cout << "- " << elements[0][1] << "c3 ";
		else if (elements[0][1] > 0)
			cout << "+ " << elements[0][1] << "c3 ";
	}
	else if (elements[0][1] != 0) {
		if (elements[0][1] == 1)
			cout << "c3 ";
		else if (elements[0][1] == -1)
			cout << "-c3 ";
		else
			cout << elements[0][1] << "c3 ";
	}
	else
		cout << 0;
	cout << "= " << answer[0] << "\n";
	for (int i = 1; i < size; ++i) {
		bool zero = true;
		bool start = false;
		if (elements[i][0] == 1) {
			cout << "c" << i + 1 << " ";
			zero = false;
		}
		else if (elements[i][0] == -1) {
			cout << "-c" << i + 1 << " ";
			zero = false;
		}
		else if (elements[i][0] != 0) {
			cout << elements[i][0] << "c" << i + 1 << " ";
			zero = false;
		}
		else
			start = true;
		for (int j = 1; j < elements[i].size(); ++j) {
			if (elements[i][j] == 1) {
				if (!start)
					cout << "+ ";
				cout << "c" << j + i + 1 << " ";
				zero = false;
			}
			else if (elements[i][j] == -1) {
				if (!start)
					cout << "- c" << j + i + 1 << " ";
				else
					cout << "-c" << j + i + 1 << " ";
				zero = false;
			}
			else if (elements[i][j] < 0) {
				if (!start)
					cout << "- " << -elements[i][j] << "c" << j + i + 1 << " ";
				else
					cout << elements[i][j] << "c" << j + i + 1 << " ";
				zero = false;
			}
			else if (elements[i][j] > 0) {
				if (!start)
					cout << "+ ";
				cout << elements[i][j] << "c" << j + i + 1 << " ";
				zero = false;
			}
			if (!zero)
				start = false;
		}
		if (zero)
			cout << 0 << " ";
		cout << "= " << answer[i] << "\n";
	}
	return;
}

void ThreeDiag::GetMatrix() {
	cout << "input size\n";
	cin >> size;
	elements.resize(size);
	for (int i = 1; i < size - 1; ++i)
		elements[i].resize(3);
	elements[0].resize(2);
	elements.back().resize(2);
	answer.resize(size);
	cout << "input matrix\n";
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < elements[i].size(); ++j)
			cin >> elements[i][j];
	return;
}

void ThreeDiag::Progonka() {
	vector<double> P(size - 1);
	P[0] = -elements[0][1] / elements[0][0];
	answer[0] /= elements[0][0];
	for (int i = 1; i < size - 1; ++i) {
		P[i] = -elements[i][2] / (elements[i][1] + elements[i][0] * P[i - 1]);
		answer[i] = (answer[i] - elements[i][0] * answer[i - 1]) /
					(elements[i][1] + elements[i][0] * P[i - 1]);
	}
	answer.back() = (answer.back() - elements.back()[0] * answer[size - 2]) /
					(elements.back()[1] + elements.back()[0] * P.back());
	for (int i = size - 2; i >= 0; --i) {
		answer[i] += P[i] * answer[i + 1];
	}
	return;
}

void ThreeDiag::CheckSolve(vector<double>& vec) {
	cout << elements[0][0] * vec[0] + elements[0][1] * vec[1] << " = " << answer[0] << "\n";
	for (int i = 1; i < size; ++i) {
		double a = 0;
		for (int j = 0; j < elements[i].size(); ++j)
			a += elements[i][j] * vec[i + j - 1];
		cout << a << " = " << answer[i] << "\n";
	}
	return;
}

vector<double>& ThreeDiag::GetAnswer() {
	return answer;
}

void ThreeDiag::PrintAnswer() {
	for (auto i : answer)
		cout << i << "\n";
	return;
}

void ThreeDiag::PrintSystem() {
	int number = 1;
	if (elements[0][0] != 0) {
		if (elements[0][0] == 1)
			cout << "x1 ";
		else if (elements[0][0] == -1)
			cout << "-x1 ";
		else
			cout << elements[0][0] << "x1 ";
		if (elements[0][1] == 1)
			cout << "+ x2 ";
		else if (elements[0][1] == -1)
			cout << "- x2 ";
		else if (elements[0][1] < 0)
			cout << "- " << elements[0][1] << "x2 ";
		else if (elements[0][1] > 0)
			cout << "+ " << elements[0][1] << "x2 ";
	}
	else if (elements[0][1] != 0) {
		if (elements[0][1] == 1)
			cout << "x2 ";
		else if (elements[0][1] == -1)
			cout << "-x2 ";
		else
			cout << elements[0][1] << "x2 ";
	}
	else
		cout << 0;
	cout << "= " << answer[0] << "\n";
	for (int i = 1; i < size; ++i) {
		bool zero = true;
		bool start = false;
		if (elements[i][0] == 1) {
			cout << "x" << i << " ";
			zero = false;
		}
		else if (elements[i][0] == -1) {
			cout << "-x" << i << " ";
			zero = false;
		}
		else if (elements[i][0] != 0) {
			cout << elements[i][0] << "x" << i << " ";
			zero = false;
		}
		else
			start = true;
		for (int j = 1; j < elements[i].size(); ++j) {
			if (elements[i][j] == 1) {
				if (!start)
					cout << "+ ";
				cout << "x" << j + i << " ";
				zero = false;
			}
			else if (elements[i][j] == -1) {
				if (!start)
					cout << "- x" << j + i << " ";
				else
					cout << "-x" << j + i << " ";
				zero = false;
			}
			else if (elements[i][j] < 0) {
				if (!start)
					cout << "- " << -elements[i][j] << "x" << j + i << " ";
				else
					cout << elements[i][j] << "x" << j + i << " ";
				zero = false;
			}
			else if (elements[i][j] > 0) {
				if (!start)
					cout << "+ ";
				cout << elements[i][j] << "x" << j + i << " ";
				zero = false;
			}
			if (!zero)
				start = false;
		}
		if (zero)
			cout << 0 << " ";
		cout << "= " << answer[i] << "\n";
	}
	return;
}

vector<double> ThreeDiag::operator [] (int index) const {
	return elements[index];
}

ostream& operator << (ostream& os, const ThreeDiag& matrix) {
	for (int i = 0; i < matrix.size; ++i) {
		for (int j = 0; j < i - 1; ++j)
			cout << "     " << 0 << " ";
		for (auto j : matrix[i])
			os << fixed << right << setfill(' ') << setw(6) << setprecision(3) << j << " ";
		for (int j = i + 2; j < matrix.size; ++j)
			cout << "     " << 0 << " ";
		os << "\n";
	}
	os << defaultfloat;
	return os;
}

int ThreeDiag::KnowSize() {
	return size;
}

void ThreeDiag::PutEndDiff(vector<vector<double>>& A) {
	size = A.size();
	elements.resize(size);
	answer.resize(size);
	for (int i = 1; i < size - 1; ++i)
		elements[i].resize(3);
	elements[0].resize(2);
	elements.back().resize(2);
	for (int i = 1; i < size - 1; ++i) {
		for (int j = 0; j < 3; ++j)
			elements[i][j] = A[j][i];
		answer[i] = A.back()[i];
	}
	elements[0][0] = A[0][0];
	elements[0][1] = A[1][0];
	elements.back()[0] = A[1].back();
	elements.back()[1] = A[2].back();
	answer[0] = A[3][0];
	answer.back() = A[3].back();
	return;
}

ThreeDiag::~ThreeDiag() {}
