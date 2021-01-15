#pragma once
#include <cmath>
#include <iostream>

class Vector
{
	int size;
	double* elem;
public:
	Vector(int dim);
	Vector(const Vector &vec);

	void addValue(unsigned int pos, double val);
	int getSize() const { return size; };
	double norm() const;

	double& operator [](int index) const { return elem[index]; };
	Vector& operator =(const Vector &vec);
	double operator *(const Vector &vec) const;
	Vector operator *(double k) const;
	Vector operator +(const  Vector &vec) const;
	Vector operator -(const  Vector &vec) const;
	Vector& operator +=(const  Vector &vec);
	Vector& operator -=(const  Vector &vec);

	friend std::ostream& operator<< (std::ostream &out, Vector &p_Vector);

	~Vector() { delete[] elem; elem = nullptr; };

};

