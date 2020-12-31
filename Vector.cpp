#include "Vector.h"

Vector::Vector(int dim) : size(dim)
{
	elem = new double[size];
	for (int i = 0; i < size; ++i)
		elem[i] = 0.;
}

Vector::Vector(const Vector &vec) : size(vec.size)
{
	elem = new double[size];
	for (int i = 0; i < size; ++i)
		elem[i] = vec[i];
}

void Vector::addValue(unsigned int pos, double val)
{
	elem[pos] += val;
}

double Vector::norm() const
{
	double tmp = 0;
	for (int i = 0; i < size; ++i)
		tmp += elem[i] * elem[i];
	return sqrt(tmp);
}

Vector& Vector::operator =(const Vector &vec)
{
	size = vec.getSize();
	for (int i = 0; i < size; ++i)
		elem[i] = vec[i];
	return (*this);
}

double Vector::operator *(const Vector &vec) const
{
	double tmp = 0;
	for (int i = 0; i < size; ++i)
		tmp += elem[i] * vec[i];
	return tmp;
}

Vector Vector::operator*(double val) const
{
	Vector tmp(*this);
	for (int i = 0; i < size; ++i)
		tmp[i] *= val;
	return tmp;
}

Vector Vector::operator +(const Vector &vec) const
{
	Vector tmp((*this).getSize());
	for (int i = 0; i < size; ++i)
		tmp[i] = elem[i] + vec[i];
	return tmp;
}

Vector Vector::operator -(const Vector &vec) const
{
	Vector tmp((*this).getSize());
	for (int i = 0; i < size; ++i)
		tmp[i] = elem[i] - vec[i];
	return tmp;
}

Vector& Vector::operator +=(const Vector &vec)
{
	for (int i = 0; i < size; ++i)
		elem[i] += vec[i];
	return (*this);
}

std::ostream& operator<< (std::ostream &out, Vector &p_Vector)
{
	for (int i = 0; i < p_Vector.size; ++i)
		out << p_Vector[i] << " ";
	out << std::endl;
	return out;
}

Vector& Vector::operator -=(const Vector &vec)
{
	for (int i = 0; i < size; ++i)
		elem[i] -= vec[i];
	return (*this);
}
