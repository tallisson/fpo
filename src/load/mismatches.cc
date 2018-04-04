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

	CalcPkB(graph);
	CalcQkB(graph);
	cout << "Mis: " << endl;
	for(int i = 0; i < m_mis.size(); i++) {
		printf("%.12f\n", m_mis(i));
	}

	return GetMis();
}

vec Mismatch::GetMis(void) {
	return m_mis;
}

}
