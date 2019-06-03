/*
 * header.h
 *
 *  Created on: Nov 20, 2017
 *      Author: Thomas
 */

#ifndef HEADER_H_
#define HEADER_H_

//PCA.c
void PCA(int k, int d, double** data);

double* multiplyMatrixVector(int r, int c1, int c2, double** WT, double* point);
double** transpose(int r, int c, double** array);
double** makeArray(int r, int c);
void freeArray(int r, double** array);
void printVector(int len, double *vector);
void printArray(int r, int c, double** array);
void setToIdentity(int r, int c, double** A);
void setAEqualB(int m, int n, double** A, double** B);
double** matrixMultiply(int r1, int c1, int r2, int c2, double** arr1, double** arr2);
double** QRalg(int m, int n, double** data);
double** makeZeroes(int r, int c);
double** outerProduct(int n, double *v1, double *v2);
void matrixAdd(int r, int c, double** arr1, double** arr2, double** arr3);
void printIVector(int len, int *vector);
double* normalize(int len, double* u);
double magnitude(int len, double* vector);


double*** decompH(int m, int n, double**data);


#endif /* HEADER_H_ */
