/*
 * QRalg.c
 *
 *  Created on: Nov 13, 2017
 *      Author: Thomas
 */

#include <stdlib.h>
#include <stdio.h>
#include "header.h"

int testConverge(int m, double* old_e, double** Anew);
int* orderE(int m, double*old_e);
void basis(int m, int* o, double** array);


double** QRalg(int m, int n, double** data) {

	double**A=data;
	//for (i=0;i<10000;i++) {
	int i=0;
	int flag=-1;
	double*old_e;
	if ((old_e=(double *)malloc(m * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	double **Q1=makeArray(m,m);
	setToIdentity(m,m,Q1);
	//while ((i<10000) && (flag !=0)) {
	while((i<535)) {
		double***QR=decompH(m,n,A);
		double**T = transpose(m,n,QR[0]);
		double**Atemp=matrixMultiply(m,n,m,n,T,A);
		double**Anew=matrixMultiply(m,n,m,n,Atemp,QR[0]);

		freeArray(m,Atemp);
		freeArray(m,A);
		A=Anew;

		double**Qttemp=matrixMultiply(m,n,m,n,Q1,QR[0]);
		setAEqualB(m,m,Q1,Qttemp);
		freeArray(m,Qttemp);
		freeArray(m,T);


		if (i!=0) {
			flag=testConverge(m, old_e, Anew);
		} else {
			int j;
			for (j=0;j<m;j++) {
				old_e[j]=A[j][j];
			}
		}

		i=i+1;
	}



	int *o=orderE(m,old_e);
	free(old_e);

	//make sure eigenvectors stored in Q1 are normalized
	double**QT=transpose(m,m,Q1);
	int a,b;
	double**ev=makeArray(m,m);
	for (a=0;a<m;a++) {
		double*col;
		if ((col=(double *)malloc(m * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
		for (b=0;b<m;b++) {
			col[b]=QT[a][b];
		}
		col=normalize(m,col);
		for (b=0;b<m;b++) {
			ev[a][b]=col[b];
		}
		free(col);
	}
	freeArray(m,QT);
	double**evf=transpose(m,m,ev);
	freeArray(m,ev);



	//order eigenvectors stored in Q1 based on eigenvalue
	basis(m,o,evf);






	freeArray(m,Q1);
	return evf;
	//return data;
}

void basis(int m, int* o, double** array) {
	double** ba = makeArray(m,m);
	int i,j;
	for (i=0;i<m;i++) {
		for (j=0;j<m;j++) {
			ba[i][j]=array[i][o[j]];
		}
	}
	setAEqualB(m,m,array,ba);
	freeArray(m,ba);
}

int* orderE(int m, double*old_e) {
	double* or_e;
	if ((or_e=(double *)malloc(m * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	int* order;
	if ((order=(int *)malloc(m * sizeof(int)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}

	int b;
	for (b=0;b<m;b++) {
		or_e[b]=old_e[b];
	}

	int i=1;
	while (i<m) {
		int j=i;
		double temp;
		//while (j>0 && old_e[j-1]>old_e[j]) {
		while (j>0 && old_e[j-1]<old_e[j]) {

			temp=old_e[j];
			old_e[j]=old_e[j-1];
			old_e[j-1]=temp;
			j=j-1;
		}
		i=i+1;
	}




	int a;
	for (a=0;a<m;a++) {
		int c,d;
		for (c=0;c<m;c++) {
			d=0;
			while (old_e[c]!=or_e[d]) {
				d=d+1;
			}
			order[c]=d;

		}
	}
	free(or_e);
	return order;
}




int testConverge(int m, double* old_e, double** Anew) {
	int flag =0;
	int i;
	for (i=0;i<m;i++) {
		if (abs(old_e[i]-Anew[i][i]) > 0.0000001) {
			flag=flag+1;
		}
		old_e[i]=Anew[i][i];
	}
	return flag;
}
