#include "CSR_Matrix.h"



#include <iostream>
#include <stdexcept>
#include <stdlib.h>

using namespace std;


CSR_Matrix::CSR_Matrix(int size) {
	rows.assign(size + 1, 0);
}

CSR_Matrix::CSR_Matrix(const CSR_Matrix& orig) {
}

CSR_Matrix::~CSR_Matrix() {
}


CSR_Matrix* CSR_Matrix::dump() {

	for (int it = 0; it < data.size(); it++)
	{
		cout << " [" << it << "] ";
		cout << data.at(it) << ", ";

	}
	cout << "\n\n";

	for (int it = 0; it < col.size(); it++)
	{
		cout << " [" << it << "] ";
		cout << col.at(it) << ", ";

	}
	cout << "\n\n";

	for (int it = 0; it < rows.size(); it++)
	{
		cout << " [" << it << "] ";
		cout << rows.at(it) << ", ";

	}
	cout << "\n\n";

	return this;
}


int CSR_Matrix::getSize() const{
	return rows.size() - 1;
}


CSR_Matrix* CSR_Matrix::setElement(int x, int y, int value) {
	if (x >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu x ");
	}
	if (y >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu y ");
	}
	int pos;
	int *positionPointer;
	positionPointer = &pos;

	if (value == 0) {
		this->removeElement(x, y);
		return this;
	}
	if (this->findPosition(x, y, positionPointer)) {
		data.at(*positionPointer) = value;
		return this;
	}
	else {
		data.insert((data.begin().operator +(*positionPointer)), value);
		col.insert((col.begin().operator +(*positionPointer)), y);

		if (x < rows.size() - 1) {
			for (int it = x + 1; it < rows.size(); it++) {
				rows.at(it)++;
			}
		}
		return this;
	}

}

CSR_Matrix* CSR_Matrix::addElement(int x, int y, int value)
{
	if (x >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu x ");
	}
	if (y >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu y ");
	}
	int pos;
	int *positionPointer;
	positionPointer = &pos;


	if (this->findPosition(x, y, positionPointer)) {
		data.at(*positionPointer) += value;
		return this;
	}
	else {
		data.insert((data.begin().operator +(*positionPointer)), value);
		col.insert((col.begin().operator +(*positionPointer)), y);

		if (x < rows.size() - 1) {
			for (int it = x + 1; it < rows.size(); it++) {
				rows.at(it)++;
			}
		}
		return this;
	}
}


int CSR_Matrix::getElement(int x, int y) {
	if (x >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu x ");
	}
	if (y >= this->getSize()) {
		throw std::invalid_argument("neplatny rozsah v argumentu y ");
	}
	int pos;
	int *positionPointer;
	positionPointer = &pos;

	if (this->findPosition(x, y, positionPointer)) {
		return data.at(*positionPointer);
	}
	else {
		return 0;
	}

}

CSR_Matrix* CSR_Matrix::removeElement(int x, int y) {
	int pos;
	int *positionPointer;
	positionPointer = &pos;

	if (this->findPosition(x, y, positionPointer)) {
		data.erase((data.begin().operator +(*positionPointer)));
		col.erase((col.begin().operator +(*positionPointer)));

		if (x < rows.size() - 1) {
			for (int it = x + 1; it < rows.size(); it++) {
				rows.at(it)--;
			}
		}
		return this;
	}
	else {
		return this;
	}
}

Vector CSR_Matrix::operator*(Vector& v) const{
	if (v.getSize() != this->getSize()) {
		throw std::invalid_argument("std::vector v musi byt rozmeru jako matice.");
	}
	Vector result = Vector(v.getSize());

	int it;
	for (it = 0; it < this->getSize(); it++) {
		for (int it2 = rows.at(it); it2 < rows.at(it + 1); it2++) {
			result[it] += data.at(it2)*v[col.at(it2)];
		}
	}

	return result;
}

CSR_Matrix* CSR_Matrix::matrixAddition(CSR_Matrix* m) {
	if (m->getSize() != this->getSize()) {
		throw std::invalid_argument("obe matice musi byt stejnych rozmeru");
	}
	for (int it = 0; it < this->getSize(); it++) {
		for (int it2 = 0; it2 < this->getSize(); it2++) {
			this->setElement(it, it2, this->getElement(it, it2) + m->getElement(it, it2));
		}
	}
	return this;
}

CSR_Matrix* CSR_Matrix::matrixMatrixProduct(CSR_Matrix* m) {
	if (m->getSize() != this->getSize()) {
		throw std::invalid_argument("obe matice musi byt stejnych rozmeru");
	}
	int sum;
	CSR_Matrix *cache;
	cache = new CSR_Matrix(this->getSize());

	for (int it = 0; it < this->getSize(); it++) {
		for (int it2 = 0; it2 < this->getSize(); it2++) {
			sum = 0;
			for (int it3 = 0; it3 < this->getSize(); it3++) {
				sum += (this->getElement(it, it3) * m->getElement(it3, it2));
			}
			cache->setElement(it, it2, sum);
		}
	}

	for (int it = 0; it < this->getSize(); it++) {
		for (int it2 = 0; it2 < this->getSize(); it2++) {
			this->setElement(it, it2, cache->getElement(it, it2));
		}
	}

	free(cache);

	return this;
}

bool CSR_Matrix::findPosition(int x, int y, int* positionPointer) {
	for (int it = rows.at(x); it < rows.at(x + 1); it++) {
		if (col.at(it) == y) {
			*positionPointer = it;
			return true;
		}
		if (col.at(it) > y) {
			*positionPointer = it;
			return false;
		}
	}
	*positionPointer = rows.at(x + 1);
	return false;

}