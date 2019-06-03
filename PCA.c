/*
 * PCA.c
 *
 *  Created on: Nov 13, 2017
 *      Author: Thomas
 */


#include <stdlib.h>
#include <stdio.h>
#include "header.h"

double* makeMeanVector(int r, int c, double** data);
double** makeCovariance(int r, int c, double* meanvector, double **data);

double var(int i, int j, int r,int c, double* meanvector,double**data);
void writeBasis(int d, double **basis);
void writePrincComp(int length, int d, double **data, double **pcaBasis);


void PCA(int k, int d, double** data) {
	//calculate mean vector
	double *meanvector=makeMeanVector(k,d,data);

	//make covariance matrix
	double **cov = makeCovariance(k,d,meanvector,data);

	//perform QRdecomposition and arrange basis into PCA basis
	double** pcaBasis=QRalg(d,d,cov);
	printf("PCA Basis\n");
	printArray(d,d,pcaBasis);
	writeBasis(d, pcaBasis);
	writePrincComp(k,d,data,pcaBasis);


	free(meanvector);
	freeArray(d,cov);
	freeArray(d,pcaBasis);
}

void writeBasis(int d, double **basis) {
	char* output ="pcaBasis.csv";
	FILE* file = fopen(output,"w");
	int i,j;
	for (i=0;i<d;i++) {
		for (j=0;j<d;j++) {
			fprintf(file,"%f,",basis[i][j]);
		}
		fprintf(file,"\n");
	}
	fclose(file);
}

void writePrincComp(int length, int d, double **data, double **pcaBasis) {
	char* output ="PC.csv";
	FILE* file = fopen(output,"w");
	double**WT=transpose(d,d,pcaBasis);
	//double**new_data=makeArray(length,d);
	int m,n,o;
	//transform data into PC space
	for (m=0;m<length;m++) {
		//get data point
		double*point;
		if ((point=(double *)malloc(d* sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
		for (n=0;n<d;n++) {
			point[n]=data[m][n];
		}
		double* new_vector= multiplyMatrixVector(d,d,d,WT,point);
		printf("\n");
		printVector(d,new_vector);
		for (o=0;o<d;o++) {
			fprintf(file,"%f,",new_vector[o]);
		}
		fprintf(file,"\n");
		free(point);

	}
	fclose(file);
	freeArray(d,WT);
}




double** makeCovariance(int r, int c, double* meanvector, double** data) {
	double** cov=makeZeroes(c,c);
	int i,j;
	for (i=0;i<c;i++) {
		for (j=0;j<c;j++) {
			cov[i][j]=var(i,j,r,c,meanvector,data);
		}
	}
	return cov;
}

double var(int i, int j, int r,int c, double* meanvector,double**data) {
	int m;
	double var=0;
	for (m=0;m<r;m++) {
		var=var+(data[m][j]-meanvector[j])*(data[m][i]-meanvector[i]);
	}
	var=var/(r-1);
	return var;
}

double* makeMeanVector(int r, int c, double** data) {
	int i,j;
	double *mean;
	if ((mean=(double *)malloc(c * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	for (i=0;i<c;i++) {
		double temp=0;
		for (j=0;j<r;j++) {
			temp=temp+data[j][i];
		}
		temp=temp/r;
		mean[i]=temp;
	}
	return mean;
}

