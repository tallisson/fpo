/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "jacobian.h"

#include <iostream>

using namespace model;

namespace load {

Jacobian::Jacobian() {

}

Jacobian::~Jacobian() {
	m_j1.clear();
	m_j2.clear();
	m_j3.clear();
	m_j4.clear();
	m_jqv.clear();
	m_matrix.clear();
}

mat Jacobian::GetMatrix() const {
	return m_matrix;
}

void Jacobian::SetMatrix(int numRows, int numCols) {
	m_matrix = zeros(numRows, numCols);
}

mat Jacobian::GetJ1() const {
	return m_j1;
}
void Jacobian::SetJ1(int numRows, int numCols) {
	m_j1 = zeros(numRows, numCols);
}

mat Jacobian::GetJ2() const {
	return m_j2;
}
void Jacobian::SetJ2(int numRows, int numCols) {
	m_j2 = zeros(numRows, numCols);
}

mat Jacobian::GetJ3() const {
	return m_j3;
}
void Jacobian::SetJ3(int numRows, int numCols) {
	m_j3 = zeros(numRows, numCols);
}

mat Jacobian::GetJ4() const {
	return m_j1;
}
void Jacobian::SetJ4(int numRows, int numCols) {
	m_j4 = zeros(numRows, numCols);
}

void Jacobian::CalcDPk(Graph* graph) {
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);
		if (busK->GetType() != Bus::SLACK) {
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int size = (int) branches.size();
			for (int j = 0; j < size; j++) {
				DBranch_t dataBranch = branches.at(j)->GetBranch();
				Bus* busM = neighbors.at(j);

				double vK, vM, aK, aM;
				vK = busK->GetVCalc();
				vM = busM->GetVCalc();
				aK = busK->GetACalc();
				aM = busM->GetACalc();

				double theta_km = aK - aM;
				double tap = dataBranch.m_tap;

				/*
				 * dPk em relação a 'ak'
				 * Jac(k-1, k-1) = -(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
				 * (s.branch.g(km)*sin(akm)-s.branch.b(km)*cos(akm)) + Jac(k-1, k-1)
				 */
				int k = busK->GetBus().m_ord;
				//if (busM->GetType () != Bus::SLACK)
				m_matrix(k, k) = -(1 / tap) * vK * vM
						* (dataBranch.m_g * sin(theta_km)
								- dataBranch.m_b * cos(theta_km))
						+ m_matrix(k, k);

				// M_J1
				m_j1(k, k) = -(1 / tap) * vK * vM
						* (dataBranch.m_g * sin(theta_km)
								- dataBranch.m_b * cos(theta_km)) + m_j1(k, k);

				/*
				 * dPk em relação a 'am' (exceto quando m for a barra slack).
				 *
				 * Jac(k-1, m-1) = (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
				 * (s.branch.g(km)*sin(akm)-s.branch.b(km)*cos(akm)) + Jac(k-1, m-1);
				 */
				if (busM->GetType() != Bus::SLACK) {
					int m = busM->GetBus().m_ord;
					m_matrix(k, m) = (1 / tap) * vK * vM
							* (dataBranch.m_g * sin(theta_km)
									- dataBranch.m_b * cos(theta_km))
							+ m_matrix(k, m);

					// M_J1
					m_j1(k, m) = (1 / tap) * vK * vM
							* (dataBranch.m_g * sin(theta_km)
									- dataBranch.m_b * cos(theta_km))
							+ m_j1(k, m);
				}

				/*
				 * dPk em relação a 'vk'.
				 * I:
				 * Jac(k-1, s.nb-1+s.bus.ordPQ(k)) =
				 * -2*s.branch.g(km)*(1/s.branch.tap(km))*s.bus.v(k) +
				 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) +
				 * Jac(k-1, s.nb-1+s.bus.ordPQ(k));
				 *
				 * II:
				 * Jac(k-1, s.nb-1+s.bus.ordPQ(k)) =
				 * -2*s.branch.g(km)*s.bus.v(k) +
				 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) +
				 * Jac(k-1, s.nb-1+s.bus.ordPQ(k));
				 */
				if (busK->GetType() == Bus::LOAD
						|| busK->GetType() == Bus::LOSS_CONTROL_REACT) {
					int index = graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ;

					if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
						m_matrix(k, index) =
								-2 * dataBranch.m_g * (1 / tap) * vK
										+ (1 / tap) * vM
												* (dataBranch.m_g
														* cos(theta_km)
														+ dataBranch.m_b
																* sin(theta_km))
										+ m_matrix(k, index);
					} else {
						m_matrix(k, index) =
								-2 * dataBranch.m_g * vK
										+ (1 / tap) * vM
												* (dataBranch.m_g
														* cos(theta_km)
														+ dataBranch.m_b
																* sin(theta_km))
										+ m_matrix(k, index);
					}
				}

				// M_J2
				int m = busK->GetBus().m_nin - 1;
				if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
					m_j2(k, m) = -2 * dataBranch.m_g * (1 / tap) * vK
							+ (1 / tap) * vM
									* (dataBranch.m_g * cos(theta_km)
											+ dataBranch.m_b * sin(theta_km))
							+ m_j2(k, m);
				} else {
					m_j2(k, m) = -2 * dataBranch.m_g * vK
							+ (1 / tap) * vM
									* (dataBranch.m_g * cos(theta_km)
											+ dataBranch.m_b * sin(theta_km))
							+ m_j2(k, m);
				}

				/*
				 * dPk em relação a 'vm'.
				 * Jac(k-1, s.nb-1+s.bus.ordPQ(m)) =
				 * (1/s.branch.tap(km))*s.bus.v(k)*
				 * (s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) + Jac(k-1, s.nb-1+s.bus.ordPQ(m));
				 */
				if (busM->GetType() == Bus::LOAD
						|| busM->GetType() == Bus::LOSS_CONTROL_REACT) {
					int index = graph->GetNumBus() - 1 + busM->GetBus().m_ordPQ;
					m_matrix(k, index) = (1 / tap) * vK
							* (dataBranch.m_g * cos(theta_km)
									+ dataBranch.m_b * sin(theta_km))
							+ m_matrix(k, index);
				}

				// M_J2
				m = busM->GetBus().m_nin - 1;
				m_j2(k, m) = (1 / tap) * vK
						* (dataBranch.m_g * cos(theta_km)
								+ dataBranch.m_b * sin(theta_km)) + m_j2(k, m);
			}
		}
	}
	std::cout.precision(11);
	std::cout.setf(ios::fixed);

	/*m_j1.raw_print(cout, "J1:");
	 std::cout << std::endl;

	 m_j2.raw_print(cout, "J2:");
	 std::cout << std::endl;*/
}

void Jacobian::CalcDQk(Graph* graph) {
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);
		if (busK->GetType() == Bus::LOAD
				|| busK->GetType() == Bus::LOSS_CONTROL_REACT) {
			int indexK = graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ;
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int size = (int) branches.size();
			for (int j = 0; j < size; j++) {
				Branch* branch = branches.at(j);
				DBranch_t dataBranch = branch->GetBranch();
				Bus* busM = neighbors.at(j);

				double vK, vM, aK, aM;
				vK = busK->GetVCalc();
				vM = busM->GetVCalc();
				aK = busK->GetACalc();
				aM = busM->GetACalc();

				double theta_km = aK - aM;
				double tap = dataBranch.m_tap;

				int k = busK->GetBus().m_ord;
				/*
				 * dQk em relaçao a 'ak'.
				 * Jac(s.nb-1+s.bus.ordPQ(k), k-1) =
				 * (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
				 * (s.branch.b(km)*sin(akm)+s.branch.g(km)*cos(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), k-1)
				 */
				if (busK->GetType() != Bus::SLACK) {
					m_matrix(indexK, k) = (1 / tap) * vK * vM
							* (dataBranch.m_b * sin(theta_km)
									+ dataBranch.m_g * cos(theta_km))
							+ m_matrix(indexK, k);
				}

				/*
				 * dQk em relaçao a 'am' (exceto quando m for a barra slack).
				 * Jac(s.nb-1+s.bus.ordPQ(k), m-1) =
				 * -(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
				 * (s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), m-1)
				 */
				if (busM->GetType() != Bus::SLACK) {
					int m = busM->GetBus().m_ord;
					m_matrix(indexK, m) = -(1 / tap) * vK * vM
							* (dataBranch.m_g * cos(theta_km)
									+ dataBranch.m_b * sin(theta_km))
							+ m_matrix(indexK, m);
				}

				/*
				 * dQk em relaçao a 'vk'
				 * I:
				 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
				 * 2*((1/s.branch.tap(km)^2)*s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
				 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
				 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
				 *
				 * II:
				 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
				 * 2*(s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
				 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
				 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k));
				 */
				/*if (busM->GetType () == Bus::LOAD ||
				 busM->GetType () == Bus::LOSS_CONTROL_REACT)*/
				if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
					m_matrix(indexK, indexK) = 2
							* (pow((1 / tap), 2) * dataBranch.m_b
									+ dataBranch.m_bsh) * vK
							- (1 / tap) * vM
									* (dataBranch.m_b * cos(theta_km)
											- dataBranch.m_g * sin(theta_km))
							+ m_matrix(indexK, indexK);
				} else {
					m_matrix(indexK, indexK) = 2
							* (dataBranch.m_b + dataBranch.m_bsh) * vK
							- (1 / tap) * vM
									* (dataBranch.m_b * cos(theta_km)
											- dataBranch.m_g * sin(theta_km))
							+ m_matrix(indexK, indexK);
				}

				/* dQk em relacao a 'vm'.
				 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m)) = -(1/s.branch.tap(km))*s.bus.v(k)*
				 * (s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m))
				 */
				if (busM->GetType() == Bus::LOAD
						|| busM->GetType() == Bus::LOSS_CONTROL_REACT) {
					int indexM = graph->GetNumBus() - 1
							+ busM->GetBus().m_ordPQ;
					m_matrix(indexK, indexM) = -(1 / tap) * vK
							* (dataBranch.m_b * cos(theta_km)
									- dataBranch.m_g * sin(theta_km))
							+ m_matrix(indexK, indexM);
				}
			}
			/*
			 * dQk em relaçao a 'vk' (continuação)
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
			 * 2*s.bus.bsh(k)*s.bus.v(k) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
			 */
			double vK;
			vK = busK->GetVCalc();
			m_matrix(indexK, indexK) = (2 * busK->GetBus().m_bsh * vK)
					+ m_matrix(indexK, indexK);
		}
	}
}

void Jacobian::CalcJDQ(Graph* graph) {
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);
		std::vector<Branch*> branches = busK->GetBranches();
		std::vector<Bus*> neighbors = busK->GetNeighbors();
		int k = busK->GetBus().m_nin - 1;
		int size = (int) branches.size();
		for (int j = 0; j < size; j++) {
			Branch* branch = branches.at(j);
			DBranch_t dataBranch = branch->GetBranch();
			Bus* busM = neighbors.at(j);

			double vK, vM, aK, aM;
			vK = busK->GetVCalc();
			vM = busM->GetVCalc();
			aK = busK->GetACalc();
			aM = busM->GetACalc();

			double theta_km = aK - aM;
			double tap = dataBranch.m_tap;

			int ordK = busK->GetBus().m_ord;
			/*
			 * dQk em relaçao a 'ak'.
			 * Jac(s.nb-1+s.bus.ordPQ(k), k-1) =
			 * (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
			 * (s.branch.b(km)*sin(akm)+s.branch.g(km)*cos(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), k-1)
			 */
			if (busK->GetType() != Bus::SLACK) {
				m_j3(k, ordK) = (1 / tap) * vK * vM
						* (dataBranch.m_b * sin(theta_km)
								+ dataBranch.m_g * cos(theta_km))
						+ m_j3(k, ordK);
			}
			/*
			 * dQk em relaçao a 'am' (exceto quando m for a barra slack).
			 * Jac(s.nb-1+s.bus.ordPQ(k), m-1) =
			 * -(1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
			 * (s.branch.g(km)*cos(akm)+s.branch.b(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), m-1)
			 */
			if (busM->GetType() != Bus::SLACK) {
				int ordM = busM->GetBus().m_ord;
				m_j3(k, ordM) = -(1 / tap) * vK * vM
						* (dataBranch.m_g * cos(theta_km)
								+ dataBranch.m_b * sin(theta_km))
						+ m_j3(k, ordM);
			}

			/*
			 * dQk em relaçao a 'vk'
			 * I:
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
			 * 2*((1/s.branch.tap(km)^2)*s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
			 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
			 *
			 * II:
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
			 * 2*(s.branch.b(km)+s.branch.bsh(km))*s.bus.v(k) -
			 * (1/s.branch.tap(km))*s.bus.v(m)*(s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) +
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k));
			 */
			if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
				m_j4(k, k) =
						2
								* (pow((1 / tap), 2) * dataBranch.m_b
										+ dataBranch.m_bsh) * vK
								- (1 / tap) * vM
										* (dataBranch.m_b * cos(theta_km)
												- dataBranch.m_g * sin(theta_km))
								+ m_j4(k, k);
			} else {
				m_j4(k, k) = 2 * (dataBranch.m_b + dataBranch.m_bsh) * vK
						- (1 / tap) * vM
								* (dataBranch.m_b * cos(theta_km)
										- dataBranch.m_g * sin(theta_km))
						+ m_j4(k, k);
			}

			/* dQk em relacao a 'vm'.
			 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m)) = -(1/s.branch.tap(km))*s.bus.v(k)*
			 * (s.branch.b(km)*cos(akm)-s.branch.g(km)*sin(akm)) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(m))
			 */
			int m = busM->GetBus().m_nin - 1;
			m_j4(k, m) = -(1 / tap) * vK
					* (dataBranch.m_b * cos(theta_km)
							- dataBranch.m_g * sin(theta_km)) + m_j4(k, m);
		}
		/*
		 * dQk em relaçao a 'vk' (continuação)
		 * Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k)) =
		 * 2*s.bus.bsh(k)*s.bus.v(k) + Jac(s.nb-1+s.bus.ordPQ(k), s.nb-1+s.bus.ordPQ(k))
		 */
		double vK;
		vK = busK->GetVCalc();
		m_j4(k, k) = (2 * busK->GetBus().m_bsh * vK) + m_j4(k, k);
	}
	std::cout.precision(11);
	std::cout.setf(ios::fixed);

	/*m_j3.raw_print(cout, "J3:");
	 std::cout << std::endl;

	 m_j4.raw_print(cout, "J4: ");
	 std::cout << std::endl;*/
}

mat Jacobian::CalcJac(Graph* graph) {
	CalcDPk(graph);
	CalcDQk(graph);
	CalcJDQ(graph);
	CalcJqv();

	return GetMatrix();
}

void Jacobian::Zeros(void) {
	int size = m_matrix.n_cols;
	m_matrix = zeros(size, size);
}

void Jacobian::Zeros(Graph* graph) {
	int nPQ = graph->GetNumPQ();
	int nPV = graph->GetNumPV();
	int numB = graph->GetNumBus();
	int num = graph->GetNumPQ() * 2 + graph->GetNumPV();
	m_matrix.clear();
	m_matrix = zeros(num, num);

	SetJ1((nPQ + nPV), (nPQ + nPV));
	SetJ2((nPQ + nPV), numB);
	SetJ3(numB, (nPQ + nPV));
	SetJ4(numB, numB);
}

void Jacobian::Zeros(int numRows, int numCols) {
	m_matrix.clear();
	m_matrix = zeros(numRows, numCols);
}

vec Jacobian::SolveSys(vec b) {
	return vec(inv(m_matrix) * -b);
}

void Jacobian::CalcJqv(void) {
	mat j1 = -(m_j1);
	mat j2 = -(m_j2);
	mat j3 = -(m_j3);
	mat j4 = -(m_j4);

	m_jqv = j4 - (j3 * inv(j1) * j2);
}

mat Jacobian::GetJqv(void) {
	return m_jqv;
}

}

