#pragma once
#include "Vector.h"
#include <vector>

class CSR_Matrix
{
	std::vector<double> data;
	std::vector<int> col;
	std::vector<int> rows;
public:
	CSR_Matrix(int size);
	CSR_Matrix(const CSR_Matrix& orig);
	virtual ~CSR_Matrix();
private:
	bool findPosition(int x, int y, int* position); 

public:
	CSR_Matrix*  dump();
	int getSize() const;
	int getElement(int x, int y);
	CSR_Matrix* setElement(int x, int y, int value);
	CSR_Matrix* addElement(int x, int y, int value);
	CSR_Matrix* removeElement(int x, int y);
	Vector operator* (Vector& v) const;
	CSR_Matrix* matrixAddition(CSR_Matrix* m);
	CSR_Matrix* matrixMatrixProduct(CSR_Matrix* m);
};

