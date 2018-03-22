/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "mismatches.h"

#include <math.h>
#include <iostream>

namespace load {

Mismatch::Mismatch() {
}

Mismatch::~Mismatch() {
	m_mis.clear();
}

/*void Mismatch::CalcPkB(Graph* graph) {
 // Balanço de potência ativa:
 for (int i = 0; i < graph->GetNumBus(); i++) {
 Bus* busK = graph->GetBus(i + 1);
 if (busK->GetType() == Bus::SLACK) {
 continue;
 }
 int k = busK->GetBus().m_ord;
 m_mis(k) = 0;

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

 //if (busK->GetBus ().m_nin == dataBranch.m_ni)
 if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
 /*mis(k-1) =
 (s.branch.g(km)*(1/s.branch.tap(km)^2)*s.bus.v(k)^2 -
 (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
 (s.branch.g(km)*cos(akm - s.branch.def(km))+s.branch.b(km)*
 sin(akm - s.branch.def(km)))) + mis(k-1);
 */
/*m_mis(k) =
 (dataBranch.m_g * 1 / pow(dataBranch.m_tap, 2)
 * pow(vK, 2)
 - (1 / dataBranch.m_tap) * vK * vM
 * (dataBranch.m_g
 * cos(
 theta_km
 - dataBranch.m_def)
 + dataBranch.m_b
 * sin(
 theta_km
 - dataBranch.m_def)))
 + m_mis(k);*/
//} else {
/*
 * mis(k-1) = (s.branch.g(km)*s.bus.v(k)^2 - (1/s.branch.tap(km))*
 * s.bus.v(k)*s.bus.v(m)*(s.branch.g(km)*cos(akm + s.branch.def(km))+
 * s.branch.b(km)*sin(akm + s.branch.def(km)))) + mis(k-1);
 */
/*std::cout << (dataBranch.m_g * pow(vK, 2) - (1 / dataBranch.m_tap) * vK * vM) << endl;
 m_mis(k) = (dataBranch.m_g * pow(vK, 2) - (1 / dataBranch.m_tap) * vK * vM
 * (dataBranch.m_g * cos(theta_km + dataBranch.m_def)
 + dataBranch.m_b * sin(theta_km + dataBranch.m_def)))
 + m_mis(k);
 }
 }

 double vK = busK->GetVCalc();
 double pgK = busK->GetPCalc();*/
/*
 * s.bus.pg(k) + s.bus.gsh(k)*s.bus.v(k)^2 - s.bus.pc(k) - mis(k-1);
 */
/*m_mis(k) = pgK + busK->GetBus().m_gsh * pow(vK, 2) - busK->GetBus().m_pc
 - m_mis(k);
 }
 }*/

void Mismatch::CalcPkB(Graph* graph) {
	// dP:
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i+1);
		if (busK->GetType() == Bus::SLACK) {
			continue;
		}
		int k = busK->GetBus().m_ord;
		m_mis(k) = busK->GetPCalc() - busK->GetBus().m_pc;
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

			m_mis(k) -= (dataBranch.m_g * pow(vK, 2) - vK * vM * (
					dataBranch.m_g * cos(theta_km) + dataBranch.m_b * sin(theta_km) ) );
		}
	}
}

/*void Mismatch::CalcQkB(Graph* graph) {
	// Balanço de potência reativa:
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);
		if (busK->GetType() == Bus::LOSS_CONTROL_REACT
				|| busK->GetType() == Bus::LOAD) {
			int index = graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ;
			m_mis(index) = 0;

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

				//if (busK->GetBus ().m_nin == dataBranch.m_ni)
				if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
					/*
					 * mis(s.nb-1+s.bus.ordPQ(k)) =
					 * (-(s.branch.b(km)*(1/(s.branch.tap(km)^2))+s.branch.bsh(km))*
					 * s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m)*
					 * (s.branch.b(km)*cos(akm - s.branch.def(km))-s.branch.g(km)*
					 * sin(akm - s.branch.def(km)))) + mis(s.nb-1+s.bus.ordPQ(k));
					 */
					/*m_mis(index) =
							(-(dataBranch.m_b * (1 / (pow(dataBranch.m_tap, 2)))
									+ dataBranch.m_bsh) * pow(vK, 2)
									+ (1 / dataBranch.m_tap) * vK * vM
											* (dataBranch.m_b
													* cos(
															theta_km
																	- dataBranch.m_def)
													- dataBranch.m_g
															* sin(
																	theta_km
																			- dataBranch.m_def)))
									+ m_mis(index);*/
				/*} else {
					/*
					 * mis(s.nb-1+s.bus.ordPQ(k)) =
					 * (-s.branch.b(km)*s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*
					 * s.bus.v(m)*(s.branch.b(km)*cos(akm + s.branch.def(km))-s.branch.g(km)*
					 * sin(akm + s.branch.def(km)))) + mis(s.nb-1+s.bus.ordPQ(k));
					 */
					/*m_mis(index) =
							(-(dataBranch.m_b + dataBranch.m_bsh) * pow(vK, 2)
									+ (1 / dataBranch.m_tap) * vK * vM
											* (dataBranch.m_b
													* cos(
															theta_km
																	+ dataBranch.m_def)
													- dataBranch.m_g
															* sin(
																	theta_km
																			+ dataBranch.m_def)))
									+ m_mis(index);
					/*m_mis(index) = (-dataBranch.m_b*pow(vK, 2) + (1/dataBranch.m_tap)*vK*
					 vM * (dataBranch.m_b * cos(theta_km + dataBranch.m_def)-dataBranch.m_g*
					 sin(theta_km + dataBranch.m_def))) + m_mis(index);*/
				//}
			//}

			/*double vK = busK->GetVCalc();
			double qgK = busK->GetQCalc();*/
			/*
			 * mis(s.nb-1+s.bus.ordPQ(k)) =
			 * s.bus.qg(k) + s.bus.bsh(k)*s.bus.v(k)^2 - s.bus.qc(k) - mis(s.nb-1+s.bus.ordPQ(k));
			 */
			/*m_mis(index) = /*busK->GetBus().m_qg*//*qgK
					+ busK->GetBus().m_bsh * pow(vK, 2) - busK->GetBus().m_qc
					- m_mis(index);
		}
	}*/
//}

void Mismatch::CalcQkB(Graph* graph) {
	// dQ:
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);

		if (busK->GetType() == Bus::LOSS_CONTROL_REACT
				|| busK->GetType() == Bus::LOAD) {
			int index = graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ;
			m_mis(index) = busK->GetQCalc() + busK->GetBus().m_bsh * pow(busK->GetVCalc(), 2) - busK->GetBus().m_qc;


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

				m_mis(index) -= ( -(dataBranch.m_b + dataBranch.m_bsh)* pow(vK, 2) + vK * vM *
								( dataBranch.m_b * cos(theta_km) - dataBranch.m_g * sin(theta_km) ) );
			}
		}
	}
}

vec Mismatch::CalcMismatches(Graph* graph) {
	m_mis = zeros<vec>(graph->GetNumPQ() * 2 + graph->GetNumPV());

	//std::cout << "CalcPkb" << std::endl;
	CalcPkB(graph);
	CalcQkB(graph);

	return GetMis();
}

vec Mismatch::GetMis(void) {
	return m_mis;
}

}
