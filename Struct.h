#pragma once
struct Node
{
	unsigned int globalID;
	double ro;
	double z;
	
	bool tetaBoundery = false;
	bool qBoundery = false;
};
struct Element
{
	unsigned int elementID;
	unsigned int nodesGlobalID[3];
	const unsigned int nodesLocalID[3] = { 0, 1, 2 };
	Node nodesCoordinates[3];
	unsigned int numberQBounderyNodes = 0;

	unsigned int getGlobalID(unsigned int localID) const { return nodesGlobalID[localID]; };
};