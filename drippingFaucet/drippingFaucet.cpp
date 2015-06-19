// drippingFaucet.cpp
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

typedef std::vector<double> DOUBLEVECTOR;

/* constants */
#define condition   0
#define rFaucet     0.5
#define Pb	        3.86
#define VInitial    0.6
#define volumeRate  1.0
#define numOfDisk   300
#define epsilon     1.0e-3

/* properties */
#define rho    1.0e3
#define sigma  0.071
#define mu	   1.0e-3
#define g      9.8
#define ds     1.0e-3
#define dt0    0.01
#define dPhi   M_PI * 0.01

/* characteristic values */
double l0;
double m0;
double P0;
double t0;
double mu0;

/* variables */
double z, z0;
double theta;
double r, rAve;
double dV, dV0;
double V;
double v0;
double dr, dr0, dr1, dr2, dr3;
double dtheta, dtheta0, dtheta1, dtheta2, dtheta3;
double dz, dz0, dz1, dz2, dz3;
int    zInitialSize, dVInitialSize, mesh, zDiskSize, dVDiskSize, rDiskSize;
int    i, j;
double x, y, phi;

/* vectors */
DOUBLEVECTOR zInitial;
DOUBLEVECTOR dVInitial;
DOUBLEVECTOR zDisk;
DOUBLEVECTOR dVDisk;
DOUBLEVECTOR rDisk;
DOUBLEVECTOR vDisk;

/* iterator */
DOUBLEVECTOR::iterator itr;

/* functions */
void RungeKuttaFehlberg(DOUBLEVECTOR &, DOUBLEVECTOR &, DOUBLEVECTOR &, DOUBLEVECTOR &);
double potentialAndViscosity(double, double, double, double, double, double, double ,double);
double potentialAndViscosityEdge(double, double, double, double ,double);
double surface(double, double, double, double, double, double);
double surfacePlus(double, double, double, double, double, double, double);
double surfaceMinus(double, double, double, double, double, double);

int _tmain(int argc, _TCHAR* argv[])
{
	if (condition != 1 && condition != 0){
		printf("error : invarid condition value. 0 or 1 is available.");
		getchar();
		return -1;
	}

	/* calculate the characteristic value */
	l0 = pow(sigma / (rho * g), 0.5);
	m0 = rho * l0 * l0 * l0;
	P0 = pow(sigma * rho * g, 0.5);
	t0 = pow(sigma / (rho * g * g * g), 0.25);
	mu0 = pow(rho * sigma * sigma * sigma / g, 0.25);

	/* set initial value */
	r = 1.0e-20;
	z = Pb;
	theta = M_PI / 2;
	dV = 0;

	zInitial.push_back(z);

	/* calculate RUnge-Kutta method */
	while (1){
		dr0 = ds * sin(theta);
		dz0 = ds * -cos(theta);
		dtheta0 = ds * (cos(theta) / r - z);

		dr1 = ds * sin(theta + dtheta0 * 0.5);
		dz1 = ds * -cos(theta + dtheta0 * 0.5);
		dtheta1 = ds * (cos(theta + dtheta0 * 0.5) / (r + dr0 * 0.5) - (z + dz0 * 0.5));

		dr2 = ds * sin(theta + dtheta1 * 0.5);
		dz2 = ds * -cos(theta + dtheta1 * 0.5);
		dtheta2 = ds * (cos(theta + dtheta1 * 0.5) / (r + dr1 * 0.5) - (z + dz1 * 0.5));

		dr3 = ds * sin(theta + dtheta2 * 0.5);
		dz3 = ds * -cos(theta + dtheta2 * 0.5);
		dtheta3 = ds * (cos(theta + dtheta2 * 0.5) / (r + dr2 * 0.5) - (z + dz2 * 0.5));

		dr = (dr0 + 2 * dr1 + 2 * dr2 + dr3) / 6;
		dz = (dz0 + 2 * dz1 + 2 * dz2 + dz3) / 6;
		dtheta = (dtheta0 + 2 * dtheta1 + 2 * dtheta2 + dtheta3) / 6;

		r += dr;
		z += dz;
		theta += dtheta;
		dV = M_PI * r * r * fabs(dz);
		V += dV;

		zInitial.push_back(z);
		dVInitial.push_back(dV);

		if (condition == 0){
			if (r > rFaucet){
				break;
			}
		}
		else if (condition == 1){
			if (V > VInitial){
				break;
			}
		}
	}

	/* discretization to disks */
	/* set z and dV initial condition */
	j = 0;
	zInitialSize = zInitial.size();
	dVInitialSize = dVInitial.size();
	mesh = (int)(zInitialSize / numOfDisk);
	for (i = 0; i < zInitialSize; i++){
		zInitial.at(i) -= zInitial.back();
	}
	for (i = 0; i < zInitialSize; i++){
		if (i % mesh == 0){
			itr = zDisk.begin();
			zDisk.insert(itr, zInitial.at(i));
			dV = 0;
			if (i + mesh < dVInitialSize){
				for (j; j < i + mesh; j++){
				dV += dVInitial.at(j);
				}
				itr = dVDisk.begin();
				dVDisk.insert(itr, dV);
			}
			else{
				for (j; j < dVInitialSize; j++){
				dV += dVInitial.at(j);
				}
				itr = dVDisk.begin();
				dVDisk.insert(itr, dV);
			}
		}
	}
	z0 = -1 * zDisk.at(0);
	dV0 = M_PI * rFaucet * rFaucet * abs(z0);
	dVDisk.at(0) += dV0;
	zDiskSize = zDisk.size();
	dVDiskSize = dVDisk.size();

	/* set r initial condition */
	for (i = 0; i < zDiskSize - 1; i++){
		rAve = pow(dVDisk.at(dVDiskSize - 1 - i) / (M_PI * abs(zDisk.at(zDiskSize - 1 - i) - zDisk.at(zDiskSize - 2 - i))), 0.5);
		itr = rDisk.begin();
		rDisk.insert(itr, rAve);
	}
	rDisk.insert(itr, rFaucet);
	rDiskSize = rDisk.size();

	/* set v initial condition */
	v0 = volumeRate / (M_PI * rFaucet * rFaucet);
	for (i = 0; i < zDiskSize; i++){
		vDisk.push_back(v0);
	}

	/* calculate Runge-Kutta-Fehlberg Method */
	RungeKuttaFehlberg(zDisk, dVDisk, rDisk, vDisk);
	for (i = 0; i < 55; i++){
		printf("%e\n", zDisk.at(i));
	}


	/* open the output file */
	FILE *fp;
	char *fname = "outputs/output.csv";
	fp = fopen(fname, "w");
	if (fp == NULL){
		printf("error : cannot open the file\n");
		return -1;
	}

	fprintf(fp, "x,y,z\n");
	for (i = 0; i < zDiskSize; i++){
		phi = 0;
		while (phi < M_PI / 2){
			x = rDisk.at(i) * cos(phi);
			y = rDisk.at(i) * sin(phi);
			fprintf(fp, "%e,%e,%e\n", x, y, zDisk.at(i));
			phi += dPhi;
		}
	}

	fclose(fp);

	/* debugger */
	if (condition == 0){
		printf("boundary condition   :radius\n");
	}
	else if (condition == 1){
		printf("boundary condition   :volume\n");
	}
	printf("radius of faucet     :%e\n", rFaucet);
	printf("initial set Volume   :%e\n", VInitial);
	printf("initial Pb           :%e\n", Pb);
	printf("r at z = z1          :%e\n", rDisk.at(1));
	printf("V at t = 0           :%e\n", V);
	printf("\n");
	printf("unit length          :%e\n", l0);
	printf("unit mass            :%e\n", m0);
	printf("unit pressure        :%e\n", P0);
	printf("unit time            :%e\n", t0);
	printf("unit viscosity       :%e\n", mu0);
	printf("\n");
	printf("Number of Disk       :%d\n", numOfDisk);
	printf("zInitialSize         :%d\n", zInitialSize);
	printf("dVInitialSize        :%d\n", dVInitialSize);
	printf("zDiskSize            :%d\n", zDiskSize);
	printf("dVDiskSize           :%d\n", dVDiskSize);
	printf("rDiskSize            :%d\n", rDiskSize);
	printf("\n");
	printf("Please press enter key.");

	getchar();

	return 0;
}


void RungeKuttaFehlberg(DOUBLEVECTOR &z, DOUBLEVECTOR &dV, DOUBLEVECTOR &r, DOUBLEVECTOR &v){
	double a2 = 1 / 4;
	double b2 = 1 / 4;
	double a3 = 3 / 8;
	double b3 = 3 / 32;
	double c3 = 9 / 32;
	double a4 = 12 / 13;
	double b4 = 1932 / 2197;
	double c4 = -7200 / 2197;
	double d4 = 7296 / 2197;
	double a5 = 1;
	double b5 = 439 / 216;
	double c5 = -8;
	double d5 = 439 / 216;
	double e5 = -845 / 4104;
	double a6 = 1 / 2;
	double b6 = -8 / 27;
	double c6 = 2;
	double d6 = -3544 / 2565;
	double e6 = 1859 / 4104;
	double f6 = -11 / 40;
	double r1 = 1 / 360;
	double r3 = -128 / 4275;
	double r4 = -2197 / 75240;
	double r5 = 1 / 50;
	double r6 = 2 / 55;
	double n1 = 25 / 216;
	double n3 = 1408 / 2565;
	double n4 = 2197 / 4014;
	double n5 = -1 / 5;

	for (i = 0; i < 10; i++){
		printf("%e\n", z.at(i));
	//	z.at(i + 1) = differential(z.at(i), dV.at(i), r.at(i), v.at(i));
	}

}

double potentialAndViscosity(double z, double zPlus, double zMinus, double dV, double dVPlus, double v, double vPlus, double vMinus){
	double potentialAndViscosity;
	double potential;
	double viscosity;
	potential = 1;
	viscosity = -3 * (-1 * ((vPlus - v) * dVPlus) / ((zPlus - z) * (zPlus - z) * dV) + (v - vMinus) / ((z - zMinus) * (z - zMinus)));
	potentialAndViscosity = potential + viscosity;
	return potentialAndViscosity;
}

double potentialAndViscosityEdge(double z, double zMinus, double v, double vMinus){
	double potentialAndViscosityEdge;
	double potential;
	double viscosity;
	potential = 1;
	viscosity = -3 * (v - vMinus) / ((z - zMinus) * (z - zMinus));
	potentialAndViscosityEdge = potential + viscosity;
	return potentialAndViscosityEdge;
}

double surface(double z, double zPlus, double zMinus, double r, double rPlus, double rMinus){
	double surface;
	surface = -M_PI * r * pow(0.25 * (zPlus - zMinus) * (zPlus - zMinus) + (r - rPlus) * (r - rPlus), 0.5) / (2 * (z - zMinus)) - M_PI * (r + rMinus) * (r - rPlus) * (r / (z - zMinus) + rPlus / (zPlus - z)) / pow((zPlus - zMinus) * (zPlus - zMinus) + (r - rPlus) * (r - rPlus), 0.5);
	return surface;
}

double surfacePlus(double z, double zPlus, double zMinus, double zMinusMinus, double r, double rPlus, double rMinus){
	double surfacePlus;
	surfacePlus = 0.5 * M_PI * (r / (z - zMinus) + rMinus / (zMinus - zMinusMinus)) * pow(0.25 * (zPlus - zMinus) * (zPlus - zMinus) + (r - rPlus) * (r - rPlus), 0.5) - M_PI * (r + rMinus) * (0.5 * (zPlus - zMinus) - (r - rPlus) * r / (z - zMinus)) / pow((zPlus - zMinus) * (zPlus - zMinus) + (r - rPlus) * (r - rPlus), 0.5);
	return surfacePlus;
}

double surfaceMinus(double z, double zPlus, double zMinus, double r, double rPlus, double rMinus){
	double surfaceMinus;
	surfaceMinus = M_PI * (r - rMinus) * (0.5 * (zPlus - zMinus) - (r - rPlus) * rPlus / (zPlus - z));
	return surfaceMinus;
}