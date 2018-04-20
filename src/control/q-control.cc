/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "q-control.h"
#include "../model/graph.h"

#include <map>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

namespace control {

QControl::QControl() {
}

QControl::~QControl() {
}

bool QControl::DoControl(Graph* graph) {
	// Cálculo da geração de potência reativa e teste de violação dos limites:
	bool update = false;

	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);

		if (busK->GetType() == Bus::SLACK
				|| busK->GetType() == Bus::GENERATION) {
			double qg = 0;
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int size = (int) branches.size();
			for (int i = 0; i < size; i++) {
				DBranch_t dataBranch = branches.at(i)->GetBranch();
				Bus* busM = neighbors.at(i);

				double vK, vM, aK, aM;
				vK = busK->GetVCalc();
				vM = busM->GetVCalc();
				aK = busK->GetACalc();
				aM = busM->GetACalc();

				double theta_km = aK - aM;

				//if (busK->GetBus ().m_nin == dataBranch.m_ni)
				if (dataBranch.m_tipo == 1 && busK->GetTap() == Bus::TAP) {
					/*
					 * s.bus.qg(k) =
					 * (-(s.branch.b(km)*(1/(s.branch.tap(km)^2))+s.branch.bsh(km))*
					 * s.bus.v(k)^2 + (1/s.branch.tap(km))*s.bus.v(k)*s.bus.v(m) *
					 * (s.branch.b(km)*cos(akm - s.branch.def(km))-s.branch.g(km)*sin(akm - s.branch.def(km)))) +
					 * s.bus.qg(k);
					 */
					double a = -(dataBranch.m_b * (1 / pow(dataBranch.m_tap, 2))
							+ dataBranch.m_bsh);
					double b = pow(vK, 2);
					double c = (1 / dataBranch.m_tap) * vK * vM;
					double d = dataBranch.m_b * cos(theta_km - dataBranch.m_def)
							- dataBranch.m_g * sin(theta_km - dataBranch.m_def);
					qg = (a * b + c * d) + qg;
				} else {
					/*
					 * s.bus.qg(k) = (-s.branch.b(km)*s.bus.v(k)^2 + (1/s.branch.tap(km)) *
					 * s.bus.v(k)*s.bus.v(m)*(s.branch.b(km)*cos(akm + s.branch.def(km)) -
					 * s.branch.g(km)*sin(akm + s.branch.def(km)))) + s.bus.qg(k);
					 */
					qg = (-(dataBranch.m_b) * pow(vK, 2) + (1 / dataBranch.m_tap) * vK * vM
						* (dataBranch.m_b * cos(theta_km + dataBranch.m_def) - dataBranch.m_g
						* sin(theta_km+ dataBranch.m_def))) + qg;
				}
			}

			double vK = busK->GetVCalc();
			// s.bus.qg(k) = - s.bus.bsh(k)*s.bus.v(k)^2 + s.bus.qc(k) + s.bus.qg(k);
			qg = -busK->GetBus().m_bsh * pow(vK, 2) + busK->GetBus().m_qc + qg;
			/*std::cout << "Bus" << busK->GetBus ().m_nin << ", Qg = " << qg << ", qgMin = " << busK->GetBus ().m_qgmin <<
			 ", qgMax = " << busK->GetBus ().m_qgmax << ", ThetaK = " << aK.Get () << std::endl;*/
			/*if (busK->GetBus().m_nin == 2) {
				std::cout << "Qg = " << qg << std::endl;
			}*/

			if (qg < busK->GetBus().m_qgmin) {
				qg = busK->GetBus().m_qgmin;
				busK->SetType(Bus::LOSS_CONTROL_REACT);
				//busK->SetVCalc(1.0);
				update = true;
			} else if (qg > busK->GetBus().m_qgmax) {;
				qg = busK->GetBus().m_qgmax;
				busK->SetType(Bus::LOSS_CONTROL_REACT);
				//busK->SetVCalc(1.0);
				update = true;
			}
			busK->SetQCalc(qg);
		}
	}
	if (update) {
		UpdateOrd(graph);
	}

	return update;
}

bool QControl::DoRestore(Graph* graph) {
	// Verificando quais barras recuperaram controle de reativo:
	bool update = false;
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* busK = graph->GetBus(i + 1);

		if (busK->GetType() == Bus::LOSS_CONTROL_REACT) {
			double vK, qgK;
			vK = busK->GetVCalc();
			qgK = busK->GetQCalc();
			if (qgK == busK->GetBus().m_qgmin && vK <= busK->GetBus().m_v) {
				//std::cout << "Min = " << vK << ", " << busK->GetBus().m_vmin << std::endl;
				busK->SetQCalc(0);
				busK->SetVCalc(busK->GetBus().m_v);
				busK->SetType(Bus::GENERATION);
				update = true;
			}

			if (qgK == busK->GetBus().m_qgmax && vK >= busK->GetBus().m_v) {
				//std::cout << "Max = " << vK << ", " << busK->GetBus().m_v << std::endl;
				busK->SetQCalc(0);
				busK->SetVCalc(busK->GetBus().m_v);
				busK->SetType(Bus::GENERATION);
				update = true;
			}
		}
	}

	if (update) {
		UpdateOrd(graph);
	}

	return update;
}

void QControl::UpdateOrd(Graph* graph) {
	int cont = 0;
	int cont1 = 0;
	int cont2 = 0;

	int npv = 0;
	int npq = 0;

	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* bus = graph->GetBus(i + 1);
		DBus_t oldBus = bus->GetBus();
		DBus_t newBus = oldBus;
		if (bus->GetType() == Bus::GENERATION) {
			newBus.m_ord = cont++;
			newBus.m_ordPV = cont1++;
			newBus.m_posPQ = i;
			bus->SetBus(newBus);
			npv++;
		}

		if (bus->GetType() == 1 || bus->GetType() == Bus::LOAD
				|| bus->GetType() == Bus::LOSS_CONTROL_REACT) {
			newBus.m_ord = cont++;
			newBus.m_ordPQ = cont2++;
			newBus.m_posPQ = i;
			bus->SetBus(newBus);
			npq++;
		}

		graph->SetNumPQ(npq);
		graph->SetNumPV(npv);
	}
}

}
