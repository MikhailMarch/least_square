//
// Created by Xiaomi on 12.02.2020.
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cstring>

#ifdef WIN32
#define GNUPLUT_NAME "C:\\Games\\gnuplot\\bin\\gnuplot -persist"
#endif 


using namespace std;
class Matrix {


public:
	int lines, colomns;

	double** a;
	Matrix T() {
		Matrix res(colomns, lines);
		for (int i = 0; i < lines; i++)
			for (int j = 0; j < colomns; j++)
				res.a[j][i] = a[i][j];
		return res;
	}

	void dispalay() {
		for (int i = 0; i < lines; i++) {
			for (int j = 0; j < colomns; j++) {
				//if (a[i][j] < 0 && a[i][j] > -0.01) a[i][j] = -a[i][j];
				if (j != colomns - 1)
					cout << setprecision(2) << fixed << a[i][j] + 0.0000001 << " ";
				else
					cout << setprecision(2) << fixed << a[i][j] + 0.0000001;
			}
			cout << "\n";
		}
	}

	Matrix(int n, int m) {
		lines = n;
		colomns = m;
		a = new double* [lines];
		for (int i = 0; i < lines; i++) {
			a[i] = new double[colomns];
		}
		for (int i = 0; i < lines; i++)
			for (int j = 0; j < colomns; j++)
				a[i][j] = 0;
	}
	Matrix(int n, int m, double** a) {
		lines = n;
		colomns = m;
		this->a = a;
	}

	void fillMatrix() {
		a = new double* [lines];
		for (int i = 0; i < lines; i++) {
			a[i] = new double[colomns];
		}
		for (int i = 0; i < lines; i++)
			for (int j = 0; j < colomns; j++)
				cin >> a[i][j];
	}
	void operator= (const Matrix& c1) {
		this->colomns = c1.colomns;
		this->lines = c1.lines;
		this->a = c1.a;
	}

	friend Matrix operator+ (Matrix const&, Matrix const&);
	friend Matrix operator- (Matrix const&, Matrix const&);
	friend Matrix operator* (Matrix const&, Matrix const&);

};

//Matrix operator= (Matrix const &c1){}
Matrix operator* (Matrix const& c1, Matrix const& c2) {
	if (c1.colomns == c2.lines) {
		Matrix res(c1.lines, c2.colomns);
		for (int i = 0; i < c1.lines; i++)
			for (int j = 0; j < c2.colomns; j++) {
				res.a[i][j] = 0;
				for (int k = 0; k < c1.colomns; k++)
					res.a[i][j] += c1.a[i][k] * c2.a[k][j];
			}
		return res;
	}
	else
		cout << "incorrect";
}

Matrix operator- (Matrix const& c1, Matrix const& c2) {
	if (c1.lines == c2.lines && c1.colomns == c2.colomns) {
		Matrix res(c1.lines, c1.colomns);
		for (int i = 0; i < c1.lines; i++)
			for (int j = 0; j < c1.colomns; j++)
				res.a[i][j] = c1.a[i][j] - c2.a[i][j];
		return res;
	}
	else {
		cout << "incorrect input matrixes" << "\n";
	}
}

Matrix operator+ (Matrix const& c1, Matrix const& c2) {
	if (c1.lines == c2.lines && c1.colomns == c2.colomns) {
		Matrix res(c1.lines, c1.colomns);
		for (int i = 0; i < c1.lines; i++)
			for (int j = 0; j < c1.colomns; j++)
				res.a[i][j] = c1.a[i][j] + c2.a[i][j];
		return res;
	}
	else {
		cout << "incorrect input matrixes" << "\n";
	}
}

class Column_vector : public Matrix {
public:
	Column_vector(int n) : Matrix(n, 1) {};
	Column_vector(int n, double** a) : Matrix(n, 1, a) {};

	void operator= (Matrix const& sm1) {
		this->colomns = sm1.colomns;
		this->lines = sm1.lines;
		this->a = sm1.a;
	}
};


class Square_matrix : public Matrix {
public:
	Square_matrix(int n, double** a) : Matrix(n, n, a) {};
	Square_matrix(int n) : Matrix(n, n) {};
	Square_matrix T() {
		Square_matrix* res = new Square_matrix(lines, a);
		Matrix* temp = res;
		*temp = temp->T();
		return *res;
	}
	void operator= (Square_matrix const& sm1) {
		this->colomns = sm1.colomns;
		this->lines = sm1.lines;
		this->a = sm1.a;
	}
	friend Square_matrix operator+ (Square_matrix const&, Square_matrix const&);
	friend Square_matrix operator- (Square_matrix const&, Square_matrix const&);
	friend Square_matrix operator* (Square_matrix const&, Square_matrix const&);
};

Square_matrix operator+ (Square_matrix const& sm1, Square_matrix const& sm2) {
	Square_matrix* res = new Square_matrix(sm1.lines, sm1.a);
	Matrix* pres = res;
	Matrix* psm2 = new Square_matrix(sm2.lines, sm2.a);
	*pres = *pres + *psm2;
	return *res;
};

Square_matrix operator- (Square_matrix const& sm1, Square_matrix const& sm2) {
	Square_matrix* res = new Square_matrix(sm1.lines, sm1.a);
	Matrix* pres = res;
	Matrix* psm2 = new Square_matrix(sm2.lines, sm2.a);
	*pres = *pres - *psm2;
	return *res;
};

Square_matrix operator* (Square_matrix const& sm1, Square_matrix const& sm2) {
	Square_matrix* res = new Square_matrix(sm1.lines, sm1.a);
	Matrix* pres = res;
	Matrix* psm2 = new Square_matrix(sm2.lines, sm2.a);
	*pres = *pres * *psm2;
	return *res;
};

class Id_matrix : public Square_matrix {
public:
	Id_matrix(int n) : Square_matrix(n) {
		lines = n;
		colomns = n;
		for (int i = 0; i < n; i++)
			a[i][i] = 1;
	}
};

class Elim_Matrix : public Id_matrix {
public:
	Elim_Matrix(int n, int el_n, int el_m, Square_matrix& sm1) : Id_matrix(n) {
		a[el_n][el_m] = -sm1.a[el_n][el_m] / sm1.a[el_m][el_m];
	}
};

class Perm_matrix : public Id_matrix {
public:
	Perm_matrix(int n, int row1, int row2, Square_matrix& sm1) : Id_matrix(n) {
		a[row1][row1] = 0;
		a[row1][row2] = 1;
		a[row2][row2] = 0;
		a[row2][row1] = 1;
	}
};

int perm(int n, int column, double** a) {
	int max = a[column][column];
	int res_row = column;
	for (int i = column + 1; i < n; i++) {
		double k = a[i][column];
		if (k < 0) k = -k;
		if (k > max) {
			max = k;
			res_row = i;
		}
	}
	return res_row;
}

Matrix Martix_concatination(Square_matrix part1, Id_matrix part2, int n) {
	Matrix res(n, 2 * n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			res.a[i][j] = part1.a[i][j];
	for (int i = 0; i < n; i++)
		for (int j = n; j < 2 * n; j++)
			res.a[i][j] = part2.a[i][j - n];
	return res;
}



Square_matrix gausse_elim(int n,  Square_matrix a) {
	Id_matrix k(n);
	bool cont = 1;
	int step = 0;
	Matrix res = Martix_concatination(a,k,n);
	Square_matrix sm = a;
	for (int i = 0; i < n - 1; i++) {
		int k = perm(n, i, sm.a);
		if (k != i) {
			Perm_matrix per(n, i, k, sm);
			res = per * res;
			sm = per * sm;
		}
		for (int j = i + 1; j < n; j++) {
			Elim_Matrix el(n, j, i, sm);
			res = el * res;
			sm = el * sm;
		}

	}
	for (int i = n - 1; i > 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			Elim_Matrix el(n, j, i, sm);
			res = el * res;
			sm = el * sm;
		}

	}
	for (int i = 0; i < n; i++) {
		Id_matrix norm(n);
		norm.a[i][i] = 1 / res.a[i][i];
		res = norm * res;
	}
	Square_matrix output(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			output.a[i][j] = res.a[i][j + n];
	return output;
}


int main() {
#ifdef WIN32
	FILE* pipe = _popen(GNUPLUT_NAME, "w");
#endif 



	
	int n = 15, deg = 3 ;
	Matrix data(n,2);
	data.fillMatrix();
	Matrix a(n, deg + 1);

	for (int i = 0; i < deg + 1; i++) {
		for (int j = 0; j < n; j++) {
			a.a[j][i] = pow(data.a[j][0], i);
		}
	}

	Matrix b(n, 1);
	for (int i = 0; i < n; i++)
		b.a[i][0] = data.a[i][1];
	//cout << "A:" << endl;
	//a.dispalay();
	Matrix temp = a.T() * a;
	//cout << "A_T*A:" << endl;
	//temp.dispalay();
	Square_matrix a_sq(deg + 1, temp.a);
	Square_matrix a_inv = gausse_elim(deg + 1, a_sq);
	//cout << "(A_T*A)^-1:" << endl;	
	//a_inv.dispalay();
	temp = a.T() * b;
	//cout << "A_T*b:" << endl;
	//temp.dispalay();
	Matrix res = a_inv * temp;
	//cout << "x~:" << endl;
	res.dispalay();
	
	if (pipe != NULL) {
		fprintf(pipe, "plot [-10:10] %f  + %f * x + %f * x**2 + %f * x**3, 'points.txt' \n", res.a[0][0], res.a[1][0], res.a[2][0], res.a[3][0]);
	}
	
}



