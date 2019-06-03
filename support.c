/*
 * support.c
 *
 *  Created on: Nov 20, 2017
 *      Author: Thomas
 */

#include <stdlib.h>
#include <stdio.h>
#include "header.h"
#include "math.h"

double* multiplyMatrixVector(int r, int c1, int c2, double** WT, double* point) {
	double* new_vector;
	if ((new_vector=(double *)malloc(r* sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}

	int i,j;
	for (i=0;i<r;i++) {
		float temp=0;
		for (j=0;j<c1;j++) {
			temp=temp+WT[i][j]*point[j];
		}
		new_vector[i]=temp;
	}
	return new_vector;

}

double* normalize(int len, double* u) {
	double mag=magnitude(len,u);
	int i;
	for (i=0;i<len;i++) {
		u[i]=u[i]/mag;
	}
	return u;
}

double magnitude(int len, double* vector) {
	int i;
	double mag = 0;
	for (i=0;i<len;i++) {
		mag = mag + vector[i]*vector[i];
	}
	mag = sqrt(mag);
	return mag;
}

void matrixAdd(int r, int c, double** arr1, double** arr2, double** arr3) {
	int i,j;
	for (i=0;i<r;i++) {
		for (j=0;j<c;j++) {
			arr3[i][j]=arr1[i][j]+arr2[i][j];
		}
	}
}

double** outerProduct(int n, double *v1, double *v2) {
	double **new_array;
	new_array=makeArray(n,n);
	int i,j;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			new_array[i][j]=v1[i]*v2[j];
		}
	}
	return new_array;
}

double** makeZeroes(int r, int c) {
	double** z=makeArray(r,c);
	int i,j;
	for (i=0;i<r;i++) {
		for (j=0;j<c;j++) {
			z[i][j]=0;
		}
	}
	return z;
}

double** transpose(int r, int c, double** array) {
	double **T = makeArray(c,r);
	int i,j;
	for (i=0;i<c;i++) {
		for (j=0;j<r;j++) {
			T[j][i]=array[i][j];
		}
	}
	return T;
}


double** matrixMultiply(int r1, int c1, int r2, int c2, double** arr1, double** arr2) {
	double **product=makeArray(c1,r2);
	if (c1!=r2) {
		printf("cannot multiply\n");
		exit(11);
	} else {
		int i,j,m;
		for (i=0;i<r1; i++) {
			for (j=0;j<c2;j++) {
				double temp_val=0;
				for (m=0;m<c1;m++) {
					temp_val=temp_val+arr1[i][m]*arr2[m][j];
				}
				product[i][j]=temp_val;
			}
		}
	}
	return product;
}

void setAEqualB(int m, int n, double** A, double** B) {
	int i,j;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			A[i][j]=B[i][j];
		}
	}
}

void setToIdentity(int r, int c, double** A) {
	int i,j;
	for (i=0;i<r;i++) {
		for (j=0;j<c;j++) {
			if (i==j) {
				A[i][j]=1;
			} else {
				A[i][j]=0;
			}
		}
	}
}

double** makeArray(int r, int c) {
	double **arr;
	int i;
	if ((arr=(double **)malloc(r * sizeof(double *)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	for (i=0;i<r;i++) {
		if ((*(arr+i)=(double *)malloc(c * sizeof(double)))==NULL) {fprintf(stderr, "malloc error\n"); exit(1);}
	}
	return arr;
}

void freeArray(int r, double** array) {
	int i;
	for (i=0;i<r;i++) {
		free(*(array+i));
	}
	free(array);
}

void printArray(int r, int c, double** array) {
	int i,j;
	for (i=0;i<r;i++) {
		for (j=0;j<c;j++) {
			printf("%.4f ",array[i][j]);
			}
		printf("\n");
	}
}

void printVector(int len, double *vector) {
	int i;
	for (i=0;i<len;i++) {
		printf("%.2f ",vector[i]);
	}
	printf("\n");
}

void printIVector(int len, int *vector) {
	int i;
	for (i=0;i<len;i++) {
		printf("%d ",vector[i]);
	}
	printf("\n");
}
