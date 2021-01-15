#include "Triangulation.h"

Triangulation Triangulation::createTriangulation(
	double ro_0, double ro_N,
	double z_0, double z_M,
	int nodeNumberRo, int nodeNumberZ)
{
	Triangulation triangulation;
	triangulation.ro_0 = ro_0;
	triangulation.ro_N = ro_N;
	triangulation.z_0 = z_0;
	triangulation.z_M = z_M;

	double roStep = (ro_N - ro_0) / (nodeNumberRo - 1);
	double zStep = (z_M - z_0) / (nodeNumberZ - 1);

	int IDCounter = 0;
	for (int roStepNumber = 0; roStepNumber < nodeNumberRo; roStepNumber++)
	{
		for (int zStepNumber = 0; zStepNumber < nodeNumberZ; zStepNumber++)
		{
			Node newNode;
			newNode.globalID = IDCounter;
			IDCounter++;
			newNode.ro = ro_0 + roStep * roStepNumber;
			newNode.z = z_0 + zStep * zStepNumber;
			if ((roStepNumber == 0) || (roStepNumber == (nodeNumberRo - 1)))
				newNode.qBoundery = true;
			if ((zStepNumber == 0) || (zStepNumber == (nodeNumberZ - 1)))
				newNode.tetaBoundery = true;
			triangulation.nodes.push_back(newNode);
		}
	}

	int elementCounter = 0;
	for (int roNodeNumber = 1; roNodeNumber < nodeNumberRo; roNodeNumber++)
	{
		for (int zNodeNumber = 0; zNodeNumber < (nodeNumberZ - 1); zNodeNumber++)
		{
			Element newElement1;

			newElement1.elementID = elementCounter;
			elementCounter++;
			newElement1.nodesGlobalID[0] = triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + zNodeNumber).globalID;
			newElement1.nodesGlobalID[1] = triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + (zNodeNumber + 1)).globalID;
			newElement1.nodesGlobalID[2] = triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + (zNodeNumber + 1)).globalID;
			newElement1.nodesCoordinates[0]
				= triangulation.nodes.at(newElement1.nodesGlobalID[0]);
			newElement1.nodesCoordinates[1]
				= triangulation.nodes.at(newElement1.nodesGlobalID[1]);
			newElement1.nodesCoordinates[2]
				= triangulation.nodes.at(newElement1.nodesGlobalID[2]);
			if (triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + zNodeNumber).qBoundery == true)
				newElement1.numberQBounderyNodes++;
			if (triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + (zNodeNumber + 1)).qBoundery == true)
				newElement1.numberQBounderyNodes++;
			if (triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + (zNodeNumber + 1)).qBoundery == true)
				newElement1.numberQBounderyNodes++;
			triangulation.elements.push_back(newElement1);
			Element newElement2;
			newElement2.elementID = elementCounter;
			elementCounter++;
			newElement2.nodesGlobalID[0] = triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + zNodeNumber).globalID;
			newElement2.nodesGlobalID[1] = triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + (zNodeNumber + 1)).globalID;
			newElement2.nodesGlobalID[2] = triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + zNodeNumber).globalID;
			newElement2.nodesCoordinates[0]
				= triangulation.nodes.at(newElement2.nodesGlobalID[0]);
			newElement2.nodesCoordinates[1]
				= triangulation.nodes.at(newElement2.nodesGlobalID[1]);
			newElement2.nodesCoordinates[2]
				= triangulation.nodes.at(newElement2.nodesGlobalID[2]);
			if (triangulation.nodes.at(
				roNodeNumber * nodeNumberZ + zNodeNumber).qBoundery == true)
				newElement2.numberQBounderyNodes++;
			if (triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + (zNodeNumber + 1)).qBoundery == true)
				newElement2.numberQBounderyNodes++;
			if (triangulation.nodes.at(
				(roNodeNumber - 1) * nodeNumberZ + zNodeNumber).qBoundery == true)
				newElement2.numberQBounderyNodes++;
			triangulation.elements.push_back(newElement2);
		}
	}
	return triangulation;
}

int Triangulation::getGlobalID(const Element& element, int localID) const
{
	return element.nodesGlobalID[localID];
}

//int Triangulation::Element::nodes_local_ID[] = { 0, 1, 2 };