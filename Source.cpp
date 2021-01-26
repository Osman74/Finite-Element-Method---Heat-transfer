#include "Matrix.h"
#include "CSR_Matrix.h"
#include "Triangulation.h"
#include <iostream>
#include <fstream>

Vector local_fk(const Element& element)
{
	Vector localVector(3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro2 = element.nodesCoordinates[1].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z2 = element.nodesCoordinates[1].z;
	double z3 = element.nodesCoordinates[2].z;

	double coef = (1 / 60)*((ro1 - ro3)*(z2 - z3) - (ro2 - ro3)*(z1 - z3));

	Matrix phi_phi_t(3, 3);
	phi_phi_t.setElement(0, 0, 3 * ro1 + 1 * ro2 + 1 * ro3);
	phi_phi_t.setElement(0, 1, 1 * ro1 + 1 * ro2 + ro3 / 2);
	phi_phi_t.setElement(0, 2, 1 * ro1 + ro2 / 2 + 1 * ro3);
	phi_phi_t.setElement(1, 0, 1 * ro1 + 1 * ro2 + ro3 / 2);
	phi_phi_t.setElement(1, 1, 1 * ro1 + 3 * ro2 + 1 * ro3);
	phi_phi_t.setElement(1, 2, ro1 / 2 + 1 * ro2 + 1 * ro3);
	phi_phi_t.setElement(2, 0, 1 * ro1 + ro2 / 2 + 1 * ro3);
	phi_phi_t.setElement(2, 1, ro1 / 2 + 1 * ro2 + 1 * ro3);
	phi_phi_t.setElement(2, 2, 1 * ro1 + 1 * ro2 + 3 * ro3);

	Vector Tk(3);
	localVector = phi_phi_t * Tk;
	localVector = localVector * coef;
	return localVector;
}

Vector local_fk_BE_1(
	const Element& element,
	const double& q_e, const double& alpha_T, const double& teta_inf)
{
	Vector localVector(3);

	double ro2 = element.nodesCoordinates[1].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z2 = element.nodesCoordinates[1].z;
	double z3 = element.nodesCoordinates[2].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 0.);
	phi_phi_t_ro.setElement(0, 1, 0.);
	phi_phi_t_ro.setElement(0, 2, 0.);
	phi_phi_t_ro.setElement(1, 0, 0.);
	phi_phi_t_ro.setElement(1, 1, 3. * ro2 + ro3);
	phi_phi_t_ro.setElement(1, 2, ro2 + ro3);
	phi_phi_t_ro.setElement(2, 0, 0.);
	phi_phi_t_ro.setElement(2, 1, ro2 + ro3);
	phi_phi_t_ro.setElement(2, 2, ro2 + 3. * ro3);

	double l23 = sqrt((ro3 - ro2) * (ro3 - ro2) + (z3 - z2) * (z3 - z2));
	double coef = (1. / 12.) * l23;

	Vector Tk(3);
	for (int i = 0; i < 3; i++)
		Tk[i] = q_e;

	localVector = phi_phi_t_ro * Tk;
	localVector = localVector * coef;

	coef = (1. / 6.)*alpha_T*teta_inf*l23;

	Vector phi_ro(3);
	phi_ro[0] = 0.;
	phi_ro[1] = 2. * ro2 + ro3;
	phi_ro[2] = ro2 + 2. * ro3;

	localVector += (phi_ro * coef);
	localVector = localVector * (-1.);
	return localVector;
}

Vector local_fk_BE_2(
	const Element& element,
	const double& q_e, const double& alpha_T, const double& teta_inf)
{
	Vector localVector(3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z3 = element.nodesCoordinates[2].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 3. * ro1 + ro3);
	phi_phi_t_ro.setElement(0, 1, 0.);
	phi_phi_t_ro.setElement(0, 2, ro1 + ro3);
	phi_phi_t_ro.setElement(1, 0, 0.);
	phi_phi_t_ro.setElement(1, 1, 0.);
	phi_phi_t_ro.setElement(1, 2, 0.);
	phi_phi_t_ro.setElement(2, 0, ro1 + ro3);
	phi_phi_t_ro.setElement(2, 1, 0.);
	phi_phi_t_ro.setElement(2, 2, ro1 + 3. * ro3);

	double l13 = sqrt((ro3 - ro1) * (ro3 - ro1) + (z3 - z1) * (z3 - z1));
	double coef = (1. / 12.) * l13;

	Vector Tk(3);
	for (int i = 0; i < 3; i++)
		Tk[i] = q_e;

	localVector = phi_phi_t_ro * Tk;
	localVector = localVector * coef;

	coef = (1. / 6.)*alpha_T*teta_inf*l13;

	Vector phi_ro(3);
	phi_ro[0] = 2. * ro1 + ro3;
	phi_ro[1] = 0.;
	phi_ro[2] = ro1 + 2. * ro3;

	localVector += (phi_ro * coef);

	localVector = localVector * (-1.);
	return localVector;
}

Vector local_fk_BE_3(
	const Element& element,
	const double& q_e, const double& alpha_T, const double& teta_inf)
{
	Vector localVector(3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro2 = element.nodesCoordinates[1].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z2 = element.nodesCoordinates[1].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 3. * ro1 + ro2);
	phi_phi_t_ro.setElement(0, 1, ro1 + ro2);
	phi_phi_t_ro.setElement(0, 2, 0.);
	phi_phi_t_ro.setElement(1, 0, ro1 + ro2);
	phi_phi_t_ro.setElement(1, 1, ro1 + 3. * ro2);
	phi_phi_t_ro.setElement(1, 2, 0.);
	phi_phi_t_ro.setElement(2, 0, 0.);
	phi_phi_t_ro.setElement(2, 1, 0.);
	phi_phi_t_ro.setElement(2, 2, 0.);

	double l12 = sqrt((ro2 - ro1) * (ro2 - ro1) + (z2 - z1) * (z2 - z1));
	double coef = (1. / 12.) * l12;

	Vector Tk(3);
	for (int i = 0; i < 3; i++)
		Tk[i] = q_e;

	localVector = phi_phi_t_ro * Tk;
	localVector = localVector * coef;

	coef = (1. / 6.)*alpha_T*teta_inf*l12;

	Vector phi_ro(3);
	phi_ro[0] = 2. * ro1 + ro2;
	phi_ro[1] = ro1 + 2. * ro2;
	phi_ro[2] = 0.;

	localVector += (phi_ro * coef);

	localVector = localVector * (-1.);
	return localVector;
}

Matrix local_Gk_BE_1(
	const Element& element, const double& alpha_T)
{
	Matrix localMatrix(3, 3);

	double ro2 = element.nodesCoordinates[1].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z2 = element.nodesCoordinates[1].z;
	double z3 = element.nodesCoordinates[2].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 0.);
	phi_phi_t_ro.setElement(0, 1, 0.);
	phi_phi_t_ro.setElement(0, 2, 0.);
	phi_phi_t_ro.setElement(1, 0, 0.);
	phi_phi_t_ro.setElement(1, 1, 3. * ro2 + ro3);
	phi_phi_t_ro.setElement(1, 2, ro2 + ro3);
	phi_phi_t_ro.setElement(2, 0, 0.);
	phi_phi_t_ro.setElement(2, 1, ro2 + ro3);
	phi_phi_t_ro.setElement(2, 2, ro2 + 3. * ro3);

	double l23 = sqrt((ro3 - ro2) * (ro3 - ro2) + (z3 - z2) * (z3 - z2));
	double coef = (1. / 12.) * alpha_T * l23 * (-1.);

	localMatrix = (phi_phi_t_ro * coef);

	return localMatrix;
}

Matrix local_Gk_BE_2(
	const Element& element, const double& alpha_T)
{
	Matrix localMatrix(3, 3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z3 = element.nodesCoordinates[2].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 3. * ro1 + ro3);
	phi_phi_t_ro.setElement(0, 1, 0.);
	phi_phi_t_ro.setElement(0, 2, ro1 + ro3);
	phi_phi_t_ro.setElement(1, 0, 0.);
	phi_phi_t_ro.setElement(1, 1, 0.);
	phi_phi_t_ro.setElement(1, 2, 0.);
	phi_phi_t_ro.setElement(2, 0, ro1 + ro3);
	phi_phi_t_ro.setElement(2, 1, 0.);
	phi_phi_t_ro.setElement(2, 2, ro1 + 3. * ro3);

	double l13 = sqrt((ro3 - ro1) * (ro3 - ro1) + (z3 - z1) * (z3 - z1));
	double coef = (1. / 12.) * alpha_T * l13 * (-1.);

	localMatrix = (phi_phi_t_ro * coef);

	return localMatrix;
}

Matrix local_Gk_BE_3(
	const Element& element, const double& alpha_T)
{
	Matrix localMatrix(3, 3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro2 = element.nodesCoordinates[1].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z2 = element.nodesCoordinates[1].z;

	Matrix phi_phi_t_ro(3, 3);
	phi_phi_t_ro.setElement(0, 0, 3. * ro1 + ro2);
	phi_phi_t_ro.setElement(0, 1, ro1 + ro2);
	phi_phi_t_ro.setElement(0, 2, 0.);
	phi_phi_t_ro.setElement(1, 0, ro1 + ro2);
	phi_phi_t_ro.setElement(1, 1, ro1 + 3. * ro2);
	phi_phi_t_ro.setElement(1, 2, 0.);
	phi_phi_t_ro.setElement(2, 0, 0.);
	phi_phi_t_ro.setElement(2, 1, 0.);
	phi_phi_t_ro.setElement(2, 2, 0.);

	double l12 = sqrt((ro2 - ro1) * (ro2 - ro1) + (z2 - z1) * (z2 - z1));
	double coef = (1. / 12.) * alpha_T * l12 * (-1.);

	localMatrix = (phi_phi_t_ro * coef);

	return localMatrix;
}

Matrix local_Gk(const Element& element, const double& lambda)
{
	Matrix localMatrix(3, 3);

	double ro1 = element.nodesCoordinates[0].ro;
	double ro2 = element.nodesCoordinates[1].ro;
	double ro3 = element.nodesCoordinates[2].ro;

	double z1 = element.nodesCoordinates[0].z;
	double z2 = element.nodesCoordinates[1].z;
	double z3 = element.nodesCoordinates[2].z;

	double vec_a[] = { ro2 - ro1, z2 - z1 };
	double vec_b[] = { ro3 - ro1, z3 - z1 };
	double vec_c[] = { ro3 - ro2, z3 - z2 };

	// Производные от функций формы по декартовым координам
	double phi1_1 = vec_a[1] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]) - vec_b[1] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
	double phi1_2 = -vec_a[0] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]) + vec_b[0] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
	double phi2_1 = vec_b[1] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
	double phi2_2 = -vec_b[0] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
	double phi3_1 = -vec_a[1] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
	double phi3_2 = vec_a[0] / (vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

	double coef = (ro1 + ro2 + ro3)*((ro1 - ro3)*(z2 - z3) - (ro2 - ro3)*(z1 - z3)) / 6;

	Matrix mLambda(2, 2);
	mLambda.setElement(0, 0, lambda);
	mLambda.setElement(1, 1, lambda);

	Matrix Bk(2, 3);
	Bk.setElement(0, 0, phi1_1);
	Bk.setElement(0, 1, phi2_1);
	Bk.setElement(0, 2, phi3_1);
	Bk.setElement(1, 0, phi1_2);
	Bk.setElement(1, 1, phi2_2);
	Bk.setElement(1, 2, phi3_2);

	localMatrix = Bk.transpose() * mLambda * Bk * coef;

	return localMatrix;
}

Vector ConjugateGradientMethod(const CSR_Matrix& globalMatrix, const Vector& globalVector)
{
	Vector x(globalVector.getSize());
	Vector r = globalVector - globalMatrix * x;
	Vector z = r;

	double epsilon = 1e-10;
	do
	{
		double alpha = (r * r) / ((globalMatrix * z) * z);
		x += z * alpha;
		Vector r_previous = r;
		r -= (globalMatrix * z) * alpha;
		double beta = (r * r) / (r_previous * r_previous);
		z = r + z * beta;
	} while (r.norm() / globalVector.norm() > epsilon);

	return x;
}

int main()
{
	double ro_0 = 300, ro_N = 320;
	double z_0 = 0, z_M = 200;
	int N = 11; //число элементов по ro
	int M = 101; //число элементов по z

	double lambda = 135.;
	double q_m = 0.;
	double ro0 = 10210.;
	double q_e = 30.;
	double alpha_T = 0;
	double teta_inf = 273.15;
	double teta_e = 298.15;

	Triangulation triangulation = Triangulation::createTriangulation(ro_0, ro_N, z_0, z_M, N, M);

	Vector globalVector(N*M);
	CSR_Matrix globalMatrix((N*M));

	for (const auto& element : triangulation.elements)
	{
		Vector localVector = local_fk(element);
		for (int i = 0; i < 3; i++)
			globalVector.addValue((element.getGlobalID(i)), localVector[i]);
	}

	for (const auto& element : triangulation.elements)
	{
		if (element.numberQBounderyNodes == 2)
		{
			Vector localVector(3);
			if (triangulation.nodes[element.nodesGlobalID[0]].qBoundery == false)
			{
				localVector = local_fk_BE_1(element, q_e, alpha_T, teta_inf);
				globalVector.addValue((element.getGlobalID(1)), localVector[1]);
				globalVector.addValue((element.getGlobalID(2)), localVector[2]);
			}
			if (triangulation.nodes[element.nodesGlobalID[1]].qBoundery == false)
			{
				localVector = local_fk_BE_2(element, q_e, alpha_T, teta_inf);
				globalVector.addValue((element.getGlobalID(0)), localVector[0]);
				globalVector.addValue((element.getGlobalID(2)), localVector[2]);
			}
			if (triangulation.nodes[element.nodesGlobalID[2]].qBoundery == false)
			{
				localVector = local_fk_BE_3(element, q_e, alpha_T, teta_inf);
				globalVector.addValue((element.getGlobalID(0)), localVector[0]);
				globalVector.addValue((element.getGlobalID(1)), localVector[1]);
			}
		}
	}

	for (auto node : triangulation.nodes)
	{
		if (node.tetaBoundery == true)
		{
			globalVector[node.globalID] = teta_e;
			for (int i = 0; i < triangulation.nodes.size(); i++)
			{
				if (i == node.globalID)
					globalMatrix.setElement(node.globalID, i, 1.);
				else
					globalMatrix.setElement(node.globalID, i, 0.);
			}
		}
	}

	for (auto element : triangulation.elements)
	{
		Matrix localMatrix(3, 3);
		localMatrix = local_Gk(element, lambda);
		for (int i = 0; i < 3; i++)
		{

			if (triangulation.nodes[element.getGlobalID(i)].tetaBoundery == false)
			{
				for (int j = 0; j < 3; j++)
				{
					if (triangulation.nodes[element.getGlobalID(j)].tetaBoundery == true)
						globalVector[element.getGlobalID(i)]
						-= localMatrix.getElements()[i][j] * globalVector[element.getGlobalID(j)];
					else
						globalMatrix.addElement(element.getGlobalID(i),
							element.getGlobalID(j), localMatrix.getElements()[i][j]);
				}
			}
		}
	}

	Vector teta = ConjugateGradientMethod(globalMatrix, globalVector);

	std::ofstream out("res.mv2");
	int nodesSize = triangulation.nodes.size();
	int elementsSize = triangulation.elements.size();
	out << nodesSize << " 3 1 teta" << std::endl;
	for (int i = 0; i < nodesSize; ++i)
	{
		out << (i + 1) << " " << triangulation.nodes[i].ro << " "
			<< triangulation.nodes[i].z << " " << 0 << " "
			<< teta[i] << std::endl;
	}

	out << elementsSize << " 3 3 BC_id mat_id mat_id_Out" << std::endl;
	for (int i = 0; i < elementsSize; ++i)
	{
		out << (i + 1) << " " << (triangulation.elements[i].nodesGlobalID[0] + 1)
			<< " " << (triangulation.elements[i].nodesGlobalID[1] + 1) << " "
			<< (triangulation.elements[i].nodesGlobalID[2] + 1) << " "
			<< 1 << " " << 1 << " " << 0 << std::endl;
	}
	
	system("pause");
	return 0;
}