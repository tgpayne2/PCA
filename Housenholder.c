/*
 * Housenholder.c
 *
 *  Created on: Nov 20, 2017
 *      Author: Thomas
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"

double ***makeQR(int r, int c, double**array);
double* copyColumn(int m, int s, double **array);
double getAlpha(int len, double *x);
double* makeE1(int len);
double* makeU(int len, double alpha, double* x, double* e1 );
double **makeQ(int r, int c, double *v);
int kronecker(int i, int j);

void updateQR(int m, int n, int s, double**Q, double***QR);


double*** decompH(int m, int n, double**data) {
	double*** QR=makeQR(m,n,data);
	int t;
	if (m-1<n) {
		t=m-1;
	} else {
		t=n;
	}
	double **A=QR[1];

	if (m<n) {
		printf("unsuitable data dimensions");
		exit(30);
	} else {
		int s;
		for (s=0;s<t;s++) {
			double* x=copyColumn(m,s,A);
			double alpha = getAlpha(m-s,x);
			double *e1=makeE1(m-s);
			double *u=makeU(m-s,alpha,x,e1);
			free(x);
			free(e1);
			double *v=normalize(m-s,u);

			double**Q=makeQ(m-s,m-s,v);
			free(v);
			//double** R=QR[1];
			updateQR(m,n,s,Q,QR);

		}

	}





	return QR;
}


void updateQR(int m, int n, int s, double**Q, double***QR) {
	double **Q_temp;
	Q_temp=makeArray(n,n);

	int i,j;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {

			if (i<s || j<s) {
				if (i==j) {
					Q_temp[i][j]=1;
				} else {
					Q_temp[i][j]=0;
				}
			}	else {
				Q_temp[i][j]=Q[i-s][j-s];
			}

		}
	}


	double**newQR1=matrixMultiply(n,n,n,n,Q_temp,QR[1]);
	setAEqualB(n,n,QR[1],newQR1);
	freeArray(n,newQR1);


	double **QT = transpose(n,n,Q_temp);
	freeArray(n,Q_temp);

	double**newQR0 = matrixMultiply(n,n,n,n,QR[0],QT);
	setAEqualB(n,n,QR[0],newQR0);
	freeArray(n,newQR0);

	freeArray(n,QT);



}
double **makeQ(int r, int c, double *v) {
	double **new_Q;
	new_Q=makeArray(r,c);

	int i,j;
	for (i=0;i<r;i++) {
		for (j=0;j<c;j++) {
			new_Q[i][j]=kronecker(i,j)-2*v[i]*v[j];
		}
	}
	return new_Q;
}


int kronecker(int i, int j) {
	int return_val;
	if (i==j) {
		return 1;
	} else {
		return 0;
	}
	return return_val;
}

double ***makeQR(int r, int c, double** data) {
	double ***arr;
	int i;
	if ((arr=(double ***)malloc(2 * sizeof(double **)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	for (i=0;i<2;i++) {
		*(arr+i)=makeArray(r,c);
	}
	setToIdentity(r,c,arr[0]);
	setAEqualB(r,c,arr[1],data);
	return arr;
}

double* copyColumn(int m, int s, double **array) {
	double* new_col;
	if ((new_col=(double *)malloc((m-s) * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	int i;
	for (i=0;i<m-s;i++) {
		new_col[i]=array[i+s][s];
	}
	return new_col;
}


double getAlpha(int len, double *x){
	int k=1;
	double mag=magnitude(len,x);
	double sign;
	if (x[k]>0) {
		sign=1;
	}else if (x[k]<0) {
		sign=-1;
	} else {
		/*
		exit(31);
		*/
		sign=0;
	}
	return sign*mag;
}

double* makeE1(int len) {
	int i;
	double *e1;
	if ((e1=(double *)malloc(len * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	for (i=0;i<len;i++) {
		if (i==0) {
			e1[i]=1;
		} else{
			e1[i]=0;
		}
	}

	return e1;
}



double* makeU(int len, double alpha, double* x, double* e1 ) {
	double *u;
	if ((u=(double *)malloc(len * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	int i;
	for (i=0;i<len;i++) {
		u[i]=x[i]-alpha*e1[i];
	}
	return u;
}


