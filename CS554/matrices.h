
#ifndef __MATRICES_H__

#define __MATRICES_H__

#include "stdafx.h"

namespace smd {

typedef struct _matrix {
	double **a;
	int rows;
	int cols;
	
	double **newarray(int rows, int cols) {
		double **newa = new double*[rows];
		for (int i = 0; i < rows; i++) {
			newa[i] = new double[cols];
		}
	}

	_matrix(int _rows, int _cols, double **_a): \
		rows(_rows), cols(_cols), a(_a) {}
	
	_matrix(int _rows, int _cols): \
		rows(_rows), cols(_cols) {
		a = newarray(rows, cols);
	}
	
	_matrix* mul(double d) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				a[i][j] *= d;
		}
		return this;
	}
	
	_matrix* div(double d) {
		return mul(1 / d);
	}
	
	_matrix* add(double d) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				a[i][j] += d;
		}
		return this;
	}
	
	_matrix* add(_matrix *b, double m) {
		if (rows == b->rows && cols == b->cols) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++)
					a[i][j] += b->a[i][j];
			}
		} else {
			printf("Matrix dimensions do not match in add()!\n");
			exit(-1);
		}
		return this;
	}
	
	_matrix* add(_matrix *b) {
		return add(b, 1);
	}

	_matrix* sub(_matrix *b) {
		return add(b, -1);
	}

	_matrix* transpose(_matrix *dest) {
		if (rows != dest->cols || cols != dest->cols) {
			printf("Matrix dimensions do not match in transpose()!\n");
			exit(-1);
		}
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				dest->a[j][i] = a[i][j];
			}
		}
		return dest;
	}
	
	_matrix* transpose() {
		double **newa = newarray(cols, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				newa[j][i] = a[i][j];
			}
		}
		delete a;
		a = newa;
		return this;
	}
	
	_matrix* mul(_matrix *b, _matrix *dest) {
		int p = b->cols;
		double **ba = b->a;
		double **ca = dest->a;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < p; j++) {
				for (int k = 0; k < cols; k++) {
					ca[i][j] += (a[i][k] * ba[k][j]);
				}
			}
		}
		return dest;
	}
	
	_matrix* setZeros() {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				a[i][j] = 0;
			}
		}
	}

} matrix;

}

#endif
