#pragma once
#include "Vector.h"
#include <iostream>

class Matrix
{
	double** elements;
	unsigned int rowsCount;
	unsigned int columnsCount;
public:
	Matrix(int n, int m);
	Matrix(const Matrix &matrix);

	double** getElements() const { return elements; };
	int getRowsCount() const { return rowsCount; };
	int getColumnsCount() const { return columnsCount; };
	void setElement(int row, int column, double element);
	void addElement(int row, int column, double element);

	Matrix transpose() const;

	Matrix& operator =(const Matrix &matrix);
	Matrix operator *(const Matrix &matrix) const;
	Vector operator *(const Vector &vector) const;
	Matrix operator *(const double &number) const;

	friend std::ostream& operator<< (std::ostream &out, Matrix &p_Matrix);

	~Matrix()
	{
		for (int i = 0; i < rowsCount; i++)
			delete[] elements[i];
		delete[] elements;
	};
};