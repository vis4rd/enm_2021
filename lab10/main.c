#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const int NX = 150;
static const int NT = 1000;
static const double DELTA = 0.1;
static const double DT = 0.05;
static const double xA = 7.5;
static const double sigma = 0.5; 
static const double xF = 2.5;

double aF(const int i, const int n, const double alpha);
void init(double* u, double* u0, double* a, const double alfa, const double beta);
void calc(const double alfa, const double beta, FILE* fp1, FILE* fp2);

int main()
{
	FILE* file00 = fopen("u00.txt", "w");
	FILE* file01 = fopen("u01.txt", "w");
	FILE* file10 = fopen("u10.txt", "w");
	FILE* filea = fopen("ua.txt", "w");
	FILE* fileE00 = fopen("e00.txt", "w");
	FILE* fileE01 = fopen("e01.txt", "w");
	FILE* fileE10 = fopen("e10.txt", "w");
	FILE* fileEa = fopen("ea.txt", "w");

	calc(0.0, 0.0, file00, fileE00);
	calc(0.0, 0.1, file01, fileE01);
	calc(0.0, 1.0, file10, fileE10);
	calc(1.0, 1.0, filea, fileEa);

	fclose(file00);
	fclose(file01);
	fclose(file10);
	fclose(filea);
	fclose(fileE00);
	fclose(fileE01);
	fclose(fileE10);
	fclose(fileEa);

	return 0;
}

double aF(const int i, const int n, const double alpha)
{
	if((i*DELTA-xF < 1e-5) && (alpha > 0.5))
	{
		return alpha * cos(50*n/NT);
	}
	else
	{
		return 0.0;
	}
}

void init(double* u, double* u0, double* a, const double alfa, const double beta)
{
	for(int i = 1; i < NX; i++)
	{
		if(alfa > 0.5)
		{
			u[i] = 0.0;
		}
		else
		{
			u[i] = exp(-pow(i*DELTA - xA, 2) / (2.0*pow(sigma, 2)));
		}
	}

	for(int i = 0; i <= NX; i++)
	{
		u0[i] = u[i];
	}

	for(int i = 1; i < NX; i++)
	{
		a[i] = (u[i+1] - (2.0*u[i]) + u[i-1]) / pow(DELTA, 2)
			 - beta * (u[i] - u0[i]) / DT
			 + aF(i, 0, alfa);
	}
}

void calc(const double alfa, const double beta, FILE* fp1, FILE* fp2)
{
	double u0[NX+1];
	double u[NX+1];
	double v[NX+1];
	double vp[NX+1];
	double a[NX+1];

	u[0] = 0.0;
	u[NX] = 0.0;
	v[0] = 0.0;
	v[NX] = 0.0;

	init(u, u0, a, alfa, beta);

	for(int n = 1; n <= NT; n++)
	{
		for(int i = 1; i < NX; i++)
		{
			vp[i] = v[i] + (0.5 * a[i] * DT);
			u0[i] = u[i];
			u[i] += (DT * vp[i]);
		}

		for(int i = 1; i < NX; i++)
		{
			a[i] = (u[i+1] - (2.0*u[i]) + u[i-1]) / pow(DELTA, 2)
				 - beta * (u[i] - u0[i]) / DT
				 + aF(i, n, alfa);
			v[i] = vp[i] + (0.5 * a[i] * DT);
		}

		double sum = 0.0;
		fprintf(fp1, "%g\n", u[0]);
		for(int i = 1; i < NX; i++)
		{
			sum += (pow(v[i], 2) + pow(0.5 * (u[i+1]-u[i-1]) / DELTA, 2));
			fprintf(fp1, "%g\n", u[i]);
		}
		fprintf(fp1, "%g\n", u[NX]);

		double E = (DELTA/4.0) * (pow( (u[1]-u[0]) / DELTA, 2)
				 + pow((u[NX]-u[NX-1]) / DELTA, 2))
				 + (0.5 * sum * DELTA);
		fprintf(fp2, "%g\n", E);
	}
}
