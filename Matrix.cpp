#include "Matrix.h"

Matrix::Matrix(int dimension1, int dimension2)
	: rowsCount(dimension1), columnsCount(dimension2)
{
	elements = new double*[rowsCount];
	for (int i = 0; i < rowsCount; i++)
	{
		elements[i] = new double[columnsCount];
		for (int j = 0; j < columnsCount; j++)
			elements[i][j] = 0.;
	}
}

Matrix::Matrix(const Matrix &matrix)
	: rowsCount(matrix.rowsCount), columnsCount(matrix.columnsCount)
{
	elements = new double*[rowsCount];
	for (int i = 0; i < rowsCount; i++)
	{
		elements[i] = new double[columnsCount];
		for (int j = 0; j < columnsCount; j++)
			elements[i][j] = matrix.getElements()[i][j];
	}
}

void Matrix::setElement(int row, int column, double element)
{
	(*this).elements[row][column] = element;
}

void Matrix::addElement(int row, int column, double element)
{
	(*this).elements[row][column] += element;
}

Matrix Matrix::transpose() const
{
	Matrix result(columnsCount, rowsCount);
	for (int i = 0; i < rowsCount; i++)
		for (int j = 0; j < columnsCount; j++)
			result.setElement(j, i, (*this).getElements()[i][j]);
	return result;
}

Matrix& Matrix::operator =(const Matrix &matrix)
{
	rowsCount = matrix.getRowsCount();
	columnsCount = matrix.getColumnsCount();
	for (int i = 0; i < rowsCount; i++)
		for (int j = 0; j < columnsCount; j++)
			elements[i][j] = matrix.getElements()[i][j];
	return (*this);
}

Matrix Matrix::operator *(const Matrix &matrix) const
{
	Matrix result(rowsCount, matrix.columnsCount);
	for (int i = 0; i < rowsCount; i++)
		for (int j = 0; j < matrix.columnsCount; j++)
			for (int k = 0; k < columnsCount; k++)
				result.addElement(
					i, j, (*this).getElements()[i][k] * matrix.getElements()[k][j]);
	return result;
}

Vector Matrix::operator *(const Vector &vector) const
{
	Vector result(rowsCount);
	for (int i = 0; i < rowsCount; i++)
		for (int j = 0; j < columnsCount; j++)
			result[i] += (*this).getElements()[i][j] * vector[j];
	return result;
}

Matrix Matrix::operator *(const double &number) const
{
	Matrix result(*this);
	for (int i = 0; i < rowsCount; i++)
		for (int j = 0; j < columnsCount; j++)
			result.setElement(i, j, (*this).getElements()[i][j] * number);
	return  result;
}

std::ostream& operator<< (std::ostream &out, Matrix &p_Matrix)
{
	for (int i = 0; i < p_Matrix.rowsCount; ++i)
	{
		for (int j = 0; j < p_Matrix.columnsCount; ++j)
			out << p_Matrix.getElements()[i][j] << " ";
		out << std::endl;
	}
	return out;
}