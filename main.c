/*
 * main.c
 *
 *  Created on: Nov 20, 2017
 *      Author: Thomas
 */
#include <stdlib.h>
#include <stdio.h>
#include "header.h"

int main() {
	//read in size of data file
	FILE *file;
	file=fopen("./Adata.csv","r");
	char line[1000];
	char *token;
	char del[]=",";
	int k=0;
	int d=0;
	while(fgets(line,sizeof(line),file)) {
		if (k==0) {
			token=strtok(line,del);
			while (token != NULL) {
				d++;
				token=strtok(NULL,del);
			}
		}
		k++;
	}
	//make data array of appropriate size and load data
	double **data=makeArray(k,d);
	rewind(file);
	k=0;
	while(fgets(line,sizeof(line),file)) {
		d=0;
		token = strtok(line,del);
		while (token != NULL) {
			data[k][d]=atof(token);
			token=strtok(NULL,del);
			d++;
		}
		k++;
	}

	PCA(k,d,data);

	free(data);
	return 0;
}
