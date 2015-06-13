// drippingFaucet.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

/* constants */
#define condition     0       /* determine B.C., 0 for radius, 1 for volume */
#define rFaucet       0.5     /* radius of the faucet */
#define Pb	          3.86	  /* pressure at the bottom of the drop */
#define VInitial      0.6     /* initial volume */
#define volumeRate    1.0     /* volume rate of faucet */ 
#define numOfDisk     300	  /* number of disk */

/* properties */
#define rho           1.0e3	       /* density of the liquid */
#define sigma         0.071        /* surface tension of the liquid */
#define mu			  1.0e-3       /* viscosity of the liquid */
#define g	          9.8          /* gravitational acceleration */
#define ds		      1.0e-3       /* minimal change of line */	 
#define dt0           0.01         /* minimal change of time */
#define dPhi          M_PI * 0.01  /* minimal change of phi for plot */

/* characteristic values */
double l0;	/* length */
double m0;	/* mass */
double P0;	/* pressure */
double t0;  /* time */
double mu0;	/* viscosity */

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

/* vectors for plotting */
vector<double> zInitial;
vector<double> dVInitial;
vector<double> zDisk;
vector<double> dVDisk;
vector<double> rDisk;
vector<double> vDisk;

/* iterator */
vector<double>::iterator itr;


int _tmain(int argc, _TCHAR* argv[])
{
	if (condition != 1 && condition != 0){
		printf("error : invarid condition value. 0 or 1 is available.\n");
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
	printf("Please press enter key.\n");

	getchar();
	/* system("output.csv"); */

	return 0;
}
