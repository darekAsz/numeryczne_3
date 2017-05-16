
#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
#define N 20
double freq[N] = { 2.160913, 2.184642, 2.208656, 2.232956, 2.257543, 2.282417, 2.307579, 2.333029,
2.358767, 2.384794, 2.411110, 2.437714, 2.464608, 2.491789, 2.519259, 2.547017, 2.575062,
2.603393, 2.632010, 2.660913 };

double s21[N] = { 0.0536239235, 0.0802905173, 0.1282381167, 0.2232705886, 0.4300353568, 0.8098952423, 0.9994597768,
0.9844939866, 0.9944789352, 0.9977955581, 0.9846761685, 0.9933531543, 0.9947537723, 0.9968944459,
0.3900713171, 0.0294907390, 0.0181270647, 0.0231763396, 0.0206478723, 0.0171024633 };

double h[N], mi[N], lambda[N], delta[N];

double macierz[N][N];
double wyniki[N];
double a[N], b[N], c[N], d[N];
double M[N];
double X[115];
double koncowa[115];
int main()
{

		////////h
	for (int j = 0; j < N - 1; j++)h[j + 1] = freq[j + 1] - freq[j];
	//////////mi
	for (int j = 1; j < N - 1; j++) mi[j] = h[j] / (h[j] + h[j + 1]);
	/////////lambda
	for (int j = 1; j < N - 1; j++)lambda[j] = h[j + 1] / (h[j] + h[j + 1]);
	////////delta
	for (int j = 1; j < N - 1; j++) delta[j] = (6 / (h[j] + h[j + 1]))*(((s21[j + 1] - s21[j]) / h[j + 1]) - ((s21[j] - s21[j - 1]) / h[j]));

	mi[0] = 0;
	lambda[0] = 0;
	delta[0] = 0;
	delta[N - 1] = 0;

	///////////macierz
	for (int i = 0; i < N; i++){
		macierz[i][i] = 2;
		wyniki[i] = delta[i];
		for (int d = i + 1; d < N; d++)		macierz[i][d] = lambda[d - 1];

		for (int s = 0; s < i; s++)macierz[i][s] = mi[s + 1];

	}

	///////zerowanie macierzy

	for (int i = N - 1; i >0; i--){
		for (int j = i - 1; j >= 0; j--){
			if (macierz[j][i] != 0){
				double mnoznik = macierz[j][i] / macierz[i][i];
				for (int d = 0; d <= i; d++)macierz[j][d] -= (macierz[i][d] * mnoznik);
				for (int y = i + 1; y < N; y++)macierz[i][y] = 0;
				wyniki[j] -= (wyniki[i] * mnoznik);
			}
		}
	}

	for (int i = 0; i < N; i++){
		for (int j = i + 1; j < N; j++){
			if (macierz[j][i] != 0){
				double mnoznik = macierz[j][i] / macierz[i][i];
				for (int d = i; d < N; d++)macierz[j][d] -= (macierz[i][d] * mnoznik);
				for (int y = 0; y < i; y++)macierz[i][y] = 0;
				wyniki[j] -= (wyniki[i] * mnoznik);
			}
		}
	}


	for (int i = 0; i < N; i++){
		for (int d = 0; d < N; d++){
			if (i != d)macierz[i][d] = 0;
		}

	}

	for (int i = 0; i < N; i++)M[i] = wyniki[i] / macierz[i][i];


	///////////wyznaczanie wspolczynnikow

	for (int j = 0; j < N - 1; j++){
		a[j] = s21[j];
		b[j] = ((s21[j + 1] - s21[j]) / h[j + 1]) - (((2 * M[j] + M[j + 1])*h[j + 1]) / 6);
		c[j] = M[j] / 2;
		d[j] = ((M[j + 1] - M[j]) / (6 * h[j + 1]));
	}



	///////tablica X[115]
	//cout << "h1: " << h[1] << "   h1/3:  " << h[1] / 3<<endl;
	int licznik = 0;
	for (int i = 0; i<115; i+=6){
		double frek = freq[licznik];
		X[i] = frek;
		X[i + 1] = frek + h[licznik+1] / 10;
		X[i + 2] = frek + h[licznik+1] / 8;
		X[i + 3] = frek + h[licznik+1] / 6;
		X[i + 4] = frek + h[licznik+1] / 4;
		X[i + 5] = frek + h[licznik+1] / 2;
		licznik++;

	}


	licznik = 0;
	for (int i = 0; i<115; i += 6){
		double frek = freq[licznik];
		koncowa[i] = a[licznik] + b[licznik] * (X[i] - freq[licznik]) + c[licznik]*pow(X[i] - freq[licznik], 2) + d[licznik] * pow(X[i] - freq[licznik], 3);
		koncowa[i+1] = a[licznik] + b[licznik] * (X[i+1] - freq[licznik]) + c[licznik] * pow(X[i+1] - freq[licznik], 2) + d[licznik] * pow(X[i+1] - freq[licznik], 3);
		koncowa[i+2] = a[licznik] + b[licznik] * (X[i+2] - freq[licznik]) + c[licznik] * pow(X[i+2] - freq[licznik], 2) + d[licznik] * pow(X[i+2] - freq[licznik], 3);
		koncowa[i+3] = a[licznik] + b[licznik] * (X[i+3] - freq[licznik]) + c[licznik] * pow(X[i+3] - freq[licznik], 2) + d[licznik] * pow(X[i+3] - freq[licznik], 3);
		koncowa[i+4] = a[licznik] + b[licznik] * (X[i+4] - freq[licznik]) + c[licznik] * pow(X[i+4] - freq[licznik], 2) + d[licznik] * pow(X[i+4] - freq[licznik], 3);
		koncowa[i+5] = a[licznik] + b[licznik] * (X[i+5] - freq[licznik]) + c[licznik] * pow(X[i+5] - freq[licznik], 2) + d[licznik] * pow(X[i+5] - freq[licznik], 3);

		licznik++;

	}
	




	//////0///
	double x = 2.182269;
	//cout << "wartosc: " << a[1] + b[1] * (x - freq[1]) + c[1] * pow(x - freq[1], 2) + d[1] * pow(x - freq[1], 3);
	double p = a[0] + b[0] * (x - freq[0]) + c[0] * pow(x - freq[0], 2) + d[0] * pow(x - freq[0], 3);
	//cout << 20 * log10(p);

	for (int i = 0; i < N; i++)cout <<i<<" :  "<< "a: " << a[i] << "   b:  " << b[i] << "  c:  " << c[i] << "  d:  " << d[i] << endl;

/*	for (int i = 0; i < N; i++){
			//cout << wyniki[i] / macierz[i][i];
			cout <<i<<" : "<< h[i];
			cout << endl;
			}*/
	
	for (int i = 0; i < 114; i++) koncowa[i] = 20 * log10(abs(koncowa[i]));


	for (int i = 0; i < 115; i++){
		//cout << setprecision(7);
		//cout << X[i] << "  ";
		//cout << setprecision(12);
		//cout << koncowa[i] << endl;

		
	}

		return 0;
}

