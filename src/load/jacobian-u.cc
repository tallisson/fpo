/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "jacobian-u.h"

namespace load {

JacobianU::JacobianU() {
	// TODO Auto-generated constructor stub
}

JacobianU::~JacobianU() {
	// TODO Auto-generated destructor stub
	m_jac.clear();
}

mat JacobianU::GetJac(void) const {
	return m_jac;
}

mat JacobianU::CalcJac(Graph* graph) {
	int numPQ = graph->GetNumPQ();
	int numPV = graph->GetNumPV();

	m_jac = zeros(numPV + 2 * numPQ, numPV + 1);

	CalcDPk(graph);
	CalcDQk(graph);

	return m_jac;
}

void JacobianU::CalcDPk(Graph* graph) {
	// Derivada de dP em relação a 'V'.
	int numBus = (int)graph->GetNumBus();
	for(int k = 0; k < numBus; k++) {
		Bus* busK = graph->GetBus(k + 1);
		if(busK->GetType() != Bus::SLACK) {
			std::vector<Branch*> branches = busK->GetBranches();
			int numBranch = (int)branches.size();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			for (int m = 0; m < numBranch; m++) {
				Bus* busM = neighbors.at(m);
				Branch* branch = branches.at(m);
				DBranch_t dataBranch = branch->GetBranch();
				double theta_km = busK->GetACalc() - busM->GetACalc();

				// dPk em relação a 'vk'.
				if(busK->GetType() == Bus::GENERATION) {
					m_jac(k-1, busK->GetOrdG()) +=
						-2*dataBranch.m_g*busK->GetVCalc() + busM->GetVCalc() *
						(dataBranch.m_g*cos(theta_km)+dataBranch.m_b*sin(theta_km));
				}

	            // dPk em relação a 'vm'.
				if(busM->GetType() != Bus::LOAD) {
					m_jac(k-1, busM->GetOrdG()) =
						busK->GetVCalc() * (dataBranch.m_g * cos(theta_km) + dataBranch.m_b * sin(theta_km));
				}
			}
		}
	}
}

void JacobianU::CalcDQk(Graph* graph) {
	// Derivada de dQ em relação a 'V'.
	int numBus = (int)graph->GetNumBus();
	for(int k = 0; k < numBus; k++) {
		Bus* busK = graph->GetBus(k + 1);
		if(busK->GetType() == Bus::LOAD) {
			std::vector<Branch*> branches = busK->GetBranches();
			int numBranch = (int)branches.size();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			for (int m = 0; m < numBranch; m++) {
				Bus* busM = neighbors.at(m);
				Branch* branch = branches.at(m);
				DBranch_t dataBranch = branch->GetBranch();
				double theta_km = busK->GetACalc() - busM->GetACalc();

	            // dQk em relaçao a 'vk'.
	            if (busK->GetType() == Bus::GENERATION) {
	            	//cout << "index = " << numBus - 1 + busK->GetBus().m_ordPV << endl;
	                m_jac(numBus - 1 + busK->GetBus().m_ordPQ, busK->GetOrdG()) +=
	                	2*(dataBranch.m_b+dataBranch.m_bsh)*busK->GetVCalc() -
						busM->GetVCalc()*(dataBranch.m_b*cos(theta_km)-dataBranch.m_g*sin(theta_km));
	            }

	            // dQk em relacao a 'vm'.
	            if (busM->GetType() != Bus::LOAD) {
	                m_jac((numBus - 1 + busK->GetBus().m_ordPQ), busM->GetOrdG()) =
	                		-busK->GetVCalc()*(dataBranch.m_b*cos(theta_km)-dataBranch.m_g*sin(theta_km));
	            }
			}

		     // dQk em relaçao a 'vk' (continuação).
		     if (busK->GetType() == Bus::GENERATION) {
		    	 m_jac (numBus - 1 + busK->GetBus().m_ordPQ, busK->GetOrdG()) +=
		    			 2*busK->GetBus().m_bsh * busK->GetVCalc();
		     }
		}
	}
}

} /* namespace load */
