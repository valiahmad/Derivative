#include <iostream>
#include <Windows.h>
#include <conio.h>
#include <math.h>
#include <stack>
using namespace std;
#define siz 100
#define col 5
#define r_ow 10
int sen[r_ow][col], x_1, x_2, L, M_m = 0, Fx = 0, ii = 0;
int targetfunc[r_ow][col], limitfunc[r_ow][col];
int dx_1[r_ow][col], dx_2[r_ow][col], dL[r_ow][col], d_2L[r_ow][col];
int row_x1 = 0, row_x2 = 0, row_l = 0, row_zdx = 0;
int matx[3][3], ans[3], D, Dx1, Dx2, Dl, zdx[r_ow][col];
int P[col], Q[col][col], length_d_2L = 0, length_zdx = 0, HB[3][3];

void reset_P_Q_HB() {//reset P,Q,HB
	int i, j;
	for (i = 0; i < r_ow; i++) {
		for (j = 0; j < col; j++) {
			HB[i][j] = Q[i][j] = P[j] = 0;
		}
	}
}
void put_HB() {//insert in HB matrix
	int i, j;
	for (i = 0; i < length_d_2L; i++) {
		HB[0][i + 1] = P[i];
		HB[i + 1][0] = P[i];
	}
	for (i = 0; i < length_zdx; i++) {
		for (j = 0; j < length_zdx; j++) {
			HB[i + 1][j + 1] = Q[i][j];
		}
	}
}
void zegondf(int i) {//second Differential of f(x)
	int j = 0;
	while (i) {
		zdx[j][1] *= zdx[j][4];
		if (zdx[j][4])zdx[j][4] -= 1;
		j++;
		i--;
	}
}
void put_Q() {//insert in Q matrix
	zegondf(row_zdx);
	int j;
	for (j = 0; j < row_zdx; j++) {
		if (zdx[j][1] != 0 && zdx[j][2] == 120) {
			if (zdx[j][0] == 45)
				Q[j][j] = -zdx[j][1];
			else Q[j][j] = zdx[j][1];
			length_zdx++;
		}
	}
}
void put_P(int i) {//insert in P matrix
	int j;
	for (j = 0; j < i; j++) {
		if (d_2L[j][1] != 0 && d_2L[j][2] == 120) {
			if (d_2L[j][0] == 45)
				P[j] = -d_2L[j][1];
			else P[j] = d_2L[j][1];
			length_d_2L++;
		}
	}
}
void set_matx() {//reset matrix
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			matx[i][j] = 0;
		}
	}
}
int solve_det(int mat[3][3]) {
	int answ;
	answ = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
		- mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
		+ mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
	return answ;
}
int det_HB() {
	int det_hb = solve_det(HB);
	return det_hb;
}
void findSolution() {
	// Matrix d using coeff as given in cramer's rule
	int d[3][3] = {
		{ matx[0][0], matx[0][1], matx[0][2] },
		{ matx[1][0], matx[1][1], matx[1][2] },
		{ matx[2][0], matx[2][1], matx[2][2] },
	};
	// Matrix d1 using coeff as given in cramer's rule
	int d1[3][3] = {
		{ ans[0], matx[0][1], matx[0][2] },
		{ ans[1], matx[1][1], matx[1][2] },
		{ ans[2], matx[2][1], matx[2][2] },
	};
	// Matrix d2 using coeff as given in cramer's rule
	int d2[3][3] = {
		{ matx[0][0], ans[0], matx[0][2] },
		{ matx[1][0], ans[1], matx[1][2] },
		{ matx[2][0], ans[2], matx[2][2] },
	};
	// Matrix d3 using coeff as given in cramer's rule
	int d3[3][3] = {
		{ matx[0][0], matx[0][1], ans[0] },
		{ matx[1][0], matx[1][1], ans[1] },
		{ matx[2][0], matx[2][1], ans[2] },
	};

	// Calculating Determinant of Matrices d, d1, d2, d3
	D = solve_det(d);
	Dx1 = solve_det(d1);
	Dx2 = solve_det(d2);
	Dl = solve_det(d3);
	// Case 1
	if (D != 0) {
		// Coeff have a unique solution. Apply Cramer's Rule
		int x1 = Dx1 / D;
		int x2 = Dx2 / D;
		int l = Dl / D; // calculating z using cramer's rule
		L = l;
		x_1 = x1;
		x_2 = x2;
	}
	// Case 2
	else {
		if (Dx1 == 0 && Dx2 == 0 && Dl == 0)
			cout << "Infinite solutions\n";
		else if (Dx1 != 0 || Dx2 != 0 || Dl != 0)
			cout << "No solutions\n";
	}
}
void put_matx() {//insert Coefficients in matrix
	int j;
	set_matx();
	for (j = 0; j < row_x1; j++) {
		if (dx_1[j][1] != 0 && dx_1[j][2] == 120 && dx_1[j][3] == 49 && dx_1[j][4] >= 1) {
			if (dx_1[j][0] == 45)
				matx[0][0] = -dx_1[j][1];
			else matx[0][0] = dx_1[j][1];
		}
		if (dx_1[j][1] != 0 && dx_1[j][2] == 120 && dx_1[j][3] == 50 && dx_1[j][4] >= 1) {
			if (dx_1[j][0] == 45)
				matx[0][1] = -dx_1[j][1];
			else matx[0][1] = dx_1[j][1];
		}
		if (dx_1[j][1] != 0 && dx_1[j][2] == 108 && dx_1[j][4] >= 1) {
			if (dx_1[j][0] == 45)
				matx[0][2] = dx_1[j][1];
			else matx[0][2] = -dx_1[j][1];
		}
		if (dx_1[j][1] != 0 && dx_1[j][4] == 0) {
			if (dx_1[j][0] == 45)
				ans[0] = dx_1[j][1];
			else ans[0] = -dx_1[j][1];
		}
	}
	for (j = 0; j < row_x2; j++) {
		if (dx_2[j][1] != 0 && dx_2[j][2] == 120 && dx_2[j][3] == 49 && dx_2[j][4] >= 1) {
			if (dx_2[j][0] == 45)
				matx[1][0] = -dx_2[j][1];
			else matx[1][0] = dx_2[j][1];
		}
		if (dx_2[j][1] != 0 && dx_2[j][2] == 120 && dx_2[j][3] == 50 && dx_2[j][4] >= 1) {
			if (dx_2[j][0] == 45)
				matx[1][1] = -dx_2[j][1];
			else matx[1][1] = dx_2[j][1];
		}
		if (dx_2[j][1] != 0 && dx_2[j][2] == 108 && dx_2[j][4] >= 1) {
			if (dx_2[j][0] == 45)
				matx[1][2] = dx_2[j][1];
			else matx[1][2] = -dx_2[j][1];
		}
		if (dx_2[j][1] != 0 && dx_2[j][4] == 0) {
			if (dx_2[j][0] == 45)
				ans[1] = dx_2[j][1];
			else ans[1] = -dx_2[j][1];
		}
	}
	for (j = 0; j < row_l; j++) {
		if (limitfunc[j][1] != 0 && limitfunc[j][2] == 120 && limitfunc[j][3] == 49 && limitfunc[j][4] >= 1) {
			if (limitfunc[j][0] == 45)
				matx[2][0] = -limitfunc[j][1];
			else matx[2][0] = limitfunc[j][1];
		}
		if (limitfunc[j][1] != 0 && limitfunc[j][2] == 120 && limitfunc[j][3] == 50 && limitfunc[j][4] >= 1) {
			if (limitfunc[j][0] == 45)
				matx[2][1] = -limitfunc[j][1];
			else matx[2][1] = limitfunc[j][1];
		}
		if (limitfunc[j][1] != 0 && limitfunc[j][2] == 108 && limitfunc[j][4] >= 1) {
			if (limitfunc[j][0] == 45)
				matx[2][2] = -limitfunc[j][1];
			else matx[2][2] = limitfunc[j][1];
		}
		if (limitfunc[j][1] != 0 && limitfunc[j][4] == 0) {
			if (limitfunc[j][0] == 45)
				ans[2] = limitfunc[j][1];
			else ans[2] = -limitfunc[j][1];
		}
	}
}
void set_sen() {
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 5; j++) {
			sen[i][j] = 48;
		}
	}
	ii = 0;
}
void cal_Fx(int i) {//calculate f(x)
	int j;
	for (j = 0; j <= i; j++) {
		if (targetfunc[j][1] != 0 && targetfunc[j][2] == 120 && targetfunc[j][3] == 49) {
			if (targetfunc[j][0] == 45)
				Fx = Fx + targetfunc[j][1] * pow(x_1, targetfunc[j][4]) * (-1);
			else Fx = Fx + targetfunc[j][1] * pow(x_1, targetfunc[j][4]);
		}
		if (targetfunc[j][1] != 0 && targetfunc[j][2] == 120 && targetfunc[j][3] == 50) {
			if (targetfunc[j][0] == 45)
				Fx = Fx + targetfunc[j][1] * pow(x_2, targetfunc[j][4]) * (-1);
			else Fx = Fx + targetfunc[j][1] * pow(x_2, targetfunc[j][4]);
		}
	}
}
void put_target(int i) {
	int j, l;
	for (l = 0; l < i; l++) {
		for (j = 0; j < 5; j++) {
			targetfunc[l][j] = sen[l][j];
		}
	}
}
void put_limit(int i) {
	int j, l;
	for (l = 0; l < i; l++) {
		for (j = 0; j < 5; j++) {
			limitfunc[l][j] = sen[l][j];
			dL[l][j] = sen[l][j];
			d_2L[l][j] = sen[l][j];
		}
	}
}
void put_dx(int i) {// differential of f(x) and L(l,x)
	int l;
	for (l = 0; l < i; l++) {
		if (sen[l][4] >= 0 && sen[l][3] == 49 && sen[l][1] != 0) {
			zdx[row_zdx][0] = dx_1[row_x1][0] = sen[l][0];
			zdx[row_zdx][1] = dx_1[row_x1][1] = sen[l][1];
			zdx[row_zdx][2] = dx_1[row_x1][2] = sen[l][2];
			zdx[row_zdx][3] = dx_1[row_x1][3] = sen[l][3];
			zdx[row_zdx][4] = dx_1[row_x1][4] = sen[l][4];
			row_x1++;
			row_zdx++;
		}
		else if (sen[l][4] >= 0 && sen[l][3] == 50 && sen[l][1] != 0) {
			zdx[row_zdx][0] = dx_2[row_x2][0] = sen[l][0];
			zdx[row_zdx][1] = dx_2[row_x2][1] = sen[l][1];
			zdx[row_zdx][2] = dx_2[row_x2][2] = sen[l][2];
			zdx[row_zdx][3] = dx_2[row_x2][3] = sen[l][3];
			zdx[row_zdx][4] = dx_2[row_x2][4] = sen[l][4];
			row_x2++;
			row_zdx++;
		}
	}
}
void put_dL(int i) {//differential of landa
	int l;
	for (l = 0; l < i; l++)dL[l][2] = 108;
	for (l = 0; l < i; l++) {
		if (dL[l][4] >= 0 && dL[l][3] == 49 && dL[l][1] != 0) {
			dx_1[row_x1][0] = dL[l][0];
			dx_1[row_x1][1] = dL[l][1];
			dx_1[row_x1][2] = dL[l][2];
			dx_1[row_x1][3] = dL[l][3];
			dx_1[row_x1][4] = dL[l][4];
			row_x1++;
		}
		else if (dL[l][4] > 0 && dL[l][3] == 50 && dL[l][1] != 0) {
			dx_2[row_x2][0] = dL[l][0];
			dx_2[row_x2][1] = dL[l][1];
			dx_2[row_x2][2] = dL[l][2];
			dx_2[row_x2][3] = dL[l][3];
			dx_2[row_x2][4] = dL[l][4];
			row_x2++;
		}
	}
}
void set_string(int i) {
	int j, l;
	for (l = 0; l < i; l++) {
		for (j = 0; j < 5; j++) {
			if (j == 4) { sen[l][j] -= 48; }
		}
	}
}
void derivative_f(int i) {
	int j = 0;
	while (i) {
		sen[j][1] *= sen[j][4];
		if (sen[j][4])sen[j][4] -= 1;
		j++;
		i--;
	}
}
void derivative_L(int i) {
	int j = 0;
	while (i) {
		d_2L[j][1] *= d_2L[j][4];
		if (d_2L[j][4])d_2L[j][4] -= 1;
		j++;
		i--;
	}
}
void sort_sen(int* st, int k) {
	stack<int>x;
	int hs, j = 0;
	while (k) {
		if (*st == 45 || *st == 43) { x.push(*st); st++; k--; }
		while (*st != 45 && *st != 43 && k) {
			x.push(*st);
			st++;
			k--;
		}
		hs = x.top();
		x.pop();
		if (x.top() == 94) { sen[ii][4] = hs; x.pop(); }
		else if ((x.top() >= 65 && x.top() <= 90) || (x.top() >= 97 && x.top() <= 122))
		{
			sen[ii][3] = hs; sen[ii][4] = 49;
		}
		else {
			hs -= 48;
			while (!x.empty() && (x.top() != 43 && x.top() != 45)) {
				hs = hs + (x.top() - 48) * (pow(10, (j + 1)));
				x.pop();
				j++;
			}
			sen[ii][1] = hs;
		}
		if (!x.empty() && (x.top() != 43 && x.top() != 45)) {
			hs = x.top();
			x.pop();
			if (hs >= 48 && hs <= 57) { sen[ii][3] = hs; }
			else { sen[ii][2] = hs; }
			if (!x.empty()) {
				hs = x.top();
				if ((hs >= 65 && hs <= 90) || (hs >= 97 && hs <= 122)) { sen[ii][2] = hs; x.pop(); }
				else {
					sen[ii][1] -= 48;
					while (!x.empty() && (x.top() != 43 && x.top() != 45)) {
						sen[ii][1] = sen[ii][1] + (x.top() - 48) * (pow(10, j));
						x.pop();
						j++;
					}
				}
				if (!x.empty() && (x.top() != 43 && x.top() != 45)) {
					sen[ii][1] -= 48;
					while (!x.empty() && (x.top() != 43 && x.top() != 45)) {
						sen[ii][1] = sen[ii][1] + (x.top() - 48) * (pow(10, j));
						x.pop();
						j++;
					}
				}
				else if (!x.empty() && sen[ii][1] == 0) { sen[ii][0] = x.top(); x.pop(); sen[ii][1] = 1; }
				if (!x.empty()) { sen[ii][0] = x.top(); x.pop(); }
				else if (x.empty() && sen[ii][0] != 45) { sen[ii][0] = 43; }
			}
			else { sen[ii][1] = 1; sen[ii][4] = 49; sen[ii][0] = 43; }
		}
		else if (!x.empty() && (x.top() == 43 || x.top() == 45)) { sen[ii][0] = x.top(); x.pop(); }
		else { sen[ii][0] = 43; }
		j = 0;
		ii++;
	}
	set_string(ii);
}
void mu_mi(int i) {//if target function was min ->multiply by minus
	if (M_m == 2) {
		for (int j = 0; j < i; j++) {
			if (sen[j][0] == 45)sen[j][0] = 43;
			else if (sen[j][0] == 43)sen[j][0] = 45;
		}
	}
}
void get_string() {
	int statement[siz], temp;
	int i;
	set_sen();
	system("cls");
	cout << "What Is The Target Function #1 Max\t#2 min\n";
	M_m = _getch() - 48;
	cout << "Enter The Target Function By This Format : Example-> -4x1^2-2x2^2+30x1+7x2\n ";
	i = 0;
	while (true) {
		temp = _getche();
		if (temp == '\r')break;
		statement[i] = temp;
		i++;
	}
	sort_sen(&statement[0], i);
	mu_mi(ii);
	put_target(ii);
	derivative_f(ii);
	put_dx(ii);
	set_sen();
	cout << "\nEnter The Limitation Function By This Format : Example-> 2x1+x2=10 -> 2x1+x2-10\n ";
	i = 0;
	while (true) {
		temp = _getche();
		if (temp == '\r')break;
		statement[i] = temp;
		i++;
	}
	sort_sen(&statement[0], i);
	row_l = ii;
	put_limit(ii);
	derivative_L(ii);
	put_dL(ii);
	put_matx();
	findSolution();
	cal_Fx(ii);
	set_sen();
	reset_P_Q_HB();
	put_P(row_l);
	put_Q();
	put_HB();
}

//////////////MAIN
int main()
{
	get_string();
	system("cls");
	int option = 1;
	cout << "Wellcome To Main Menu\nChoose The Option\n";
	cout << "#1 Print x_1 , x_2 , Landa , f(x)\t#2 Print P\t#3 Print Q\t#4 Print HB\n";
	cout << "#5 Print Det.HB(Result)\t\t#6 Exit\n";
	while (option != '6') {
		option = _getch();
		if (option == 49) {
			cout << "\n(x_1 = " << x_1 << ")(x_2 = " << x_2 << ")(Landa = " << L << ")(f(x) = " << Fx << ")\n";
		}
		else if (option == 50) {
			cout << "\nP = [";
			for (int i = 0; i < length_d_2L; i++)
				cout << " " << P[i]; cout << " ]\n";
		}
		else if (option == 51) {
			cout << "\nQ = [";
			for (int i = 0; i < length_zdx; i++) {
				cout << endl;
				for (int j = 0; j < length_zdx; j++) {
					cout << Q[i][j] << " ";
				}
			}
			cout << " ]\n";
		}
		else if (option == 52) {
			cout << "\nHB = [";
			for (int i = 0; i < length_zdx + 1; i++) {
				cout << endl;
				for (int j = 0; j < length_zdx + 1; j++) {
					cout << HB[i][j] << " ";
				}
			}
			cout << " ]\n";
		}
		else if (option == 53) {
			if (det_HB() >= 0) {
				cout << "\n(Det.HB = " << det_HB() << ") It Is Max";
			}
			else if (det_HB() < 0) {
				cout << "\n(Det.HB = " << det_HB() << ") It Is min";
			}
		}
	}


	/*THIS IS FOR TEST VARIABLE
	cout << "\nthis is zdx ";
	for (int i = 0; i < 10; i++) {
		cout << endl;
		for (int j = 0; j < 5; j++) {
			if (j == 1 || j == 4) { cout << " " << zdx[i][j]; continue; }
			cout << " " << (char)zdx[i][j];
		}
	}
	cout << "\nthis is d_2L  ";
	for (int i = 0; i < 10; i++) {
		cout << endl;
		for (int j = 0; j < 5; j++) {
			if (j == 1 || j == 4) { cout << " " << d_2L[i][j]; continue; }
			cout << " " << (char)d_2L[i][j];
		}
	}
	cout << "\nthis is dx 1 ";
	for (int i = 0; i < 10; i++) {
		cout << endl;
		for (int j = 0; j < 5; j++) {
			if (j == 1||j==4) { cout << " " << dx_1[i][j]; continue; }
			cout << " "<< (char)dx_1[i][j];
		}
	}
	cout << "\nthis is dx  2";
	for (int i = 0; i < 10; i++) {
		cout << endl;
		for (int j = 0; j < 5; j++) {
			if (j == 1 || j == 4) { cout << " " << dx_2[i][j]; continue; }
			cout << " " << (char)dx_2[i][j];
		}
	}
	cout << "\nthis is limit function ";
	for (int i = 0; i < 10; i++) {
		cout << endl;
		for (int j = 0; j < 5; j++) {
			if (j == 1 || j == 4) { cout << " " << limitfunc[i][j]; continue; }
			cout << " " << (char)limitfunc[i][j];
		}
	}
	for (int i = 0; i < 3; i++) {
		cout << endl;
		for (int j = 0; j < 3; j++) {
			cout << " " << matx[i][j];
		}
		cout << " = " << ans[i];
	}
	cout << "*************length of d_2L**************";
	for (int i = 0; i < length_d_2L; i++)
		cout << P[i];*/

	return 0;
}