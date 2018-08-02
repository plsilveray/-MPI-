#include "10个函数的计算.h"
int NVARS = 15;
int dim = 30;
double f1(double x[])
{
	double rt1 = 0;
	int i;
	for (i = 0; i<dim; i++) {
		rt1 += (x[i] * x[i]);
	}
	return rt1;
}

double f2(double x[])
{
	double rt1, rt2;
	int i;
	rt1 = 0, rt2 = 1;
	for (i = 0; i<dim; i++) {
		if (x[i]>0) {
			rt1 += x[i];
			rt2 *= x[i];
		}
		else {
			rt1 -= x[i];
			rt2 *= -x[i];
		}
	}
	rt1 += rt2;
	return rt1;
}

double f3(double x[])
{
	double rt1, rt2;
	int i, j;
	rt1 = 0;
	for (i = 0; i<dim; i++) {
		rt2 = 0;
		for (j = 0; j <= i; j++) {
			rt2 += x[j];
		}
		rt1 += (rt2*rt2);
	}
	return rt1;
}

double f4(double x[])
{
	double rt1;
	int i;

	rt1 = 0;
	for (i = 0; i<dim; i++) {
		if (x[i]>0 && x[i] > rt1) {
			rt1 = x[i];
		}
		if (x[i]<0 && -x[i] > rt1) {
			rt1 = -x[i];
		}
	}
	return rt1;
}

double f5(double x[])
{
	double rt1, rt2;
	int i;
	rt1 = 0, rt2 = 0;
	for (i = 0; i<dim - 1; i++) {
		rt1 += 100 * ((x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i])) + (x[i] - 1)*(x[i] - 1);
	}
	return rt1;
}

double f6(double x[])
{
	double rt1;
	int i;
	rt1 = 0;
	for (i = 0; i<dim; i++) {
		rt1 += (x[i] + 0.5)*(x[i] + 0.5);
	}
	return rt1;
}

double f7(double x[])
{
	double rt1;
	int i;
	rt1 = 0;
	for (i = 0; i<dim; i++) {
		rt1 += (i + 1)*x[i] * x[i] * x[i] * x[i];
	}
	rt1 += rand()*1.0 / RAND_MAX;
	return rt1;
}

double f8(double x[])
{
	double rt1, rt2;
	int i;
	rt1 = 0;
	for (i = 0; i<dim; i++) {
		if (x[i]<0)rt2 = -x[i];
		else rt2 = x[i];
		rt1 += -x[i] * sin(sqrt(rt2));
	}
	return rt1;
}

double f9(double x[])
{
	double rt1;
	int i;
	double pi;
	rt1 = 0;
	pi = acos(-1.0);
	for (i = 0; i<dim; i++) {
		rt1 += x[i] * x[i] - 10 * cos(2 * pi*x[i]) + 10;
	}
	return rt1;
}

double f10(double x[])
{
	double rt1, rt2;
	int i;
	double pi;
	rt1 = 0;
	rt2 = 0;
	pi = acos(-1.0);
	for (i = 0; i<dim; i++) {
		rt1 += x[i] * x[i];
		rt2 += cos(2 * pi*x[i]);
	}
	rt1 = -20 * exp(-0.2*sqrt(rt1 / dim)) - exp(rt2 / dim) + 20 + exp(1);
	return rt1;
}