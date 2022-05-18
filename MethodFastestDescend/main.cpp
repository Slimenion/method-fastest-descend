#include <iostream>
using namespace std;
//f(x) =x^2 + 7y^2 - xy + x  [1,1;1,1]

double f(double* x) {
	return pow(x[0], 2) + 7 * pow(x[1], 2) - x[0] * x[1] + x[2] / 2 * pow(x[0] + x[1] - 2, 2);
}

double ftx(double* x, double t, double* grad, int sizeX) {
	double* xt = new double[sizeX];
	for (int i = 0; i < sizeX; i++) {
		xt[i] = x[i] - t * grad[i];
	}
	return f(xt);
}

double* findGradMassiv(double* x, double numOfValues, double eps = 0.000001) {
	double* ans = new double[numOfValues];
	for (int i = 0; i < numOfValues; i++) {
		double* xWithEps = new double[numOfValues];
		for (int i = 0; i < numOfValues; i++) {
			xWithEps[i] = x[i];
		}
		xWithEps[i] += eps;
		ans[i] = (f(xWithEps) - f(x)) / eps;
	}
	return ans;
}

double findNorm(double* grad,double numOfValues) {
	double sumSquare = 0;
	for (int i = 0; i < numOfValues; i++) {
		sumSquare += pow(grad[i], 2);
	}
	sumSquare = pow(sumSquare, 0.5);
	return pow(sumSquare, 0.5);
}

void printAnswer(double* x,double numOfValues, string reasonForStop) {
	cout << "Программа завершила работу из-за:" << endl << reasonForStop << endl;
	cout << "Функция достигает минимума в точке:" << endl;
	cout << "x(";
	for (int i = 0; i < numOfValues; i++) {
		if (i != numOfValues - 1) {
			cout << x[i] << ";";
		}
		else {
			cout << x[i] << ")" << endl;
		}
	}
	cout << "Функция имеет значение в точке:" << endl;
	cout << f(x) << endl;
}

double* findIntervalMethodSvenna(double* x, double numOfValues, double t, double h) {
	double a, b;
	//Метод Свенна
	double xt = t;
	a = xt - h;//берем точку [0;0]
	b = xt + h;
	double fxPlusT = ftx(x,b,findGradMassiv(x,numOfValues),numOfValues), fxMinusT = ftx(x, a, findGradMassiv(x, numOfValues), numOfValues), fx = ftx(x, xt, findGradMassiv(x, numOfValues), numOfValues);
	if (fxMinusT <= fx && fx >= fxPlusT) {
		cout << "Функция не унимодальна" << endl;
	}
	else {
		if (fxMinusT >= fx && fx <= fxPlusT) {
			//cout << "a0 = " << a << " ; b0 = " << b << endl;
		}
		else {
			if (fxMinusT >= fx && fx >= fxPlusT) {
				a = xt;
				xt = xt + h;
			}
			else {
				h = -h;
				b = xt;
				xt = xt + h;
			}
			int k = 1;
			double x1;
			while (true) {
				x1 = xt + pow(2, k) * h;
				double deltafx = ftx(x, x1, findGradMassiv(x, numOfValues), numOfValues) - ftx(x, xt, findGradMassiv(x, numOfValues), numOfValues);
				if (deltafx < 0) {
					if (h > 0) {
						a = xt;
					}
					if (h < 0) {
						b = xt;
					}
					xt = x1;
					k++;
				}
				else {
					if (h > 0) {
						b = x1;
					}
					if (h < 0) {
						a = x1;
					}
					break;
				}

			}
		}
	}
	cout << "Итоговый промежуток [ " << a << "; " << b << " ]" << endl;
	double ans[2] = {a, b};
	return ans;
}

double findMinZolotoeSech(double* x, double numOfValues, double* interval) {
	double ak, zolotoeSceh = 0.382, bk, delta = 0.2, eps = 0.5, xmin;
	ak = interval[0];
	bk = interval[1];
	int k = 0;
	double yk, zk;
	double l = 2 * eps;
	double value_fYK, value_fZK;
	yk = ak + zolotoeSceh * (bk - ak);
	zk = ak + bk - yk;
	value_fYK = ftx(x, yk, findGradMassiv(x, numOfValues), numOfValues);
	value_fZK = ftx(x, zk, findGradMassiv(x, numOfValues), numOfValues);
	while (true)
	{

		if (value_fYK <= value_fZK) {
			bk = zk;
			zk = yk;
			yk = ak + bk - yk;
			k++;
			value_fZK = value_fYK;
			value_fYK = ftx(x, yk, findGradMassiv(x, numOfValues), numOfValues);

		}
		else {
			ak = yk;
			yk = zk;
			zk = ak + bk - zk;
			k++;
			value_fYK = value_fZK;
			value_fZK = ftx(x, zk, findGradMassiv(x, numOfValues), numOfValues);
		}
		if (abs(bk - ak) <= l) {
			xmin = (ak + bk) / 2;
			break;
		}

	}
	return xmin;
}

int main()
{
	setlocale(LC_ALL, "ru");
	int numOfValues;
	double eps1,eps2, M;
	cout << "Введите количество переменных в функции" << endl;
	cin >> numOfValues;
	double* x = new double[numOfValues];
	double* x1 = new double[numOfValues];
	double* deltaX = new double[numOfValues];
	double t = 0;
	for (int i = 0; i < numOfValues; i++) {
		cout << "Введите x" << i << endl;
		cin >> x[i];
	}
	cout << "Введите eps1" << endl;
	cin >> eps1;
	cout << "Введите eps2" << endl;
	cin >> eps2;
	cout << "Введите макисмальное число итераций" << endl;
	cin >> M;
	int k = 0;
	bool flagK = false;
	while (true) {
		double* grad = findGradMassiv(x, numOfValues);
		cout << "Градиент в точке равен:" << endl;
		for (int i = 0; i < numOfValues; i++) {
			cout << grad[i] << " ";
		}
		cout << endl;
		cout << "Норма градиентов равна:" << endl << findNorm(grad, numOfValues) << endl;
		if (findNorm(grad, numOfValues) < eps1) {
			printAnswer(x, numOfValues, "Норма градиентов меньше eps1");
			break;
		}
		if (k >= M) {
			printAnswer(x, numOfValues, "Превышено количество максимальных итераций");
			break;
		}
		else {
			t = findMinZolotoeSech(x, numOfValues, findIntervalMethodSvenna(x, numOfValues,t, 0.1));
			cout <<"t: "<< t << endl;
			cout << "x1: " << endl;
			for (int i = 0; i < numOfValues; i++) {
				x1[i] = x[i] - t * grad[i];
				cout << x1[i] << " ";
				deltaX[i] = x1[i] - x[i];
			}
			cout << endl;
			cout << "Норма разницы: " << findNorm(deltaX, numOfValues) << " Модуль значений функции: " << abs(f(x1) - f(x)) << endl;
			if (findNorm(deltaX, numOfValues) < eps2 && abs(f(x1) - f(x)) < eps2) {
				if (flagK) {
					x = x1;
					printAnswer(x, numOfValues, "Норма разницы векторов и разница значений функций дважды меньше, чем поставленные оценки");
					break;
				}
				flagK = true;
			}
			else {
				flagK = false;
			}
			for (int i = 0; i < numOfValues; i++) {
				x[i] = x1[i];
			}
			cout << "Конец " << k << " шага" << endl<<endl<<endl;
			k++;
		}

	}

}