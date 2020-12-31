#pragma once
#include <vector>
#include "Struct.h"

class Triangulation
{	
public:
	Triangulation() {};
	~Triangulation() {};
	static Triangulation createTriangulation(double ro_0, double ro_N,
		double z_0, double z_M,
		int nodeNumberRo, int nodeNumberZ);

	int getGlobalID(const Element& element, int localID) const;

	double ro_0;
	double ro_N;
	double z_0;
	double z_M;

	std::vector<Node> nodes;
	std::vector<Element> elements;
};

