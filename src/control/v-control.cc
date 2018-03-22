/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "v-control.h"

namespace control {

VControl::VControl() {
}

VControl::~VControl() {
}

Bus* VControl::MaxDsv(Graph* graph) {
	Bus* maxBus = graph->GetBus(1);
	double maxDsv = maxBus->GetDsv();
	for (int i = 1; i <= graph->GetNumBus(); i++) {
		Bus* bus = graph->GetBus(i);
		if (bus->GetType() != Bus::LOAD
				&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			continue;
		}

		if (bus->GetDsv() > maxDsv) {
			maxBus = bus;
			maxDsv = bus->GetDsv();
		}

	}

	return maxBus;
}

int VControl::MaxV(Graph* graph, arma::vec vt, Bus* modBus) {
	double maxQ = 0;
	int maxBus = 0;
	int size = (int)vt.n_elem;
	for (int i = 0; i < size; ++i) {
		Bus* crtBus = graph->GetBus(i + 1);
		if (crtBus->GetType() == Bus::LOAD
				|| crtBus->GetType() == Bus::LOSS_CONTROL_REACT) {
			continue;
		}

		crtBus->SetCrt(vt(i));
		double v;
		v = crtBus->GetVCalc();

		if (v == Bus::MAX_VOLTAGE_ONS
				&& modBus->GetStatus() == Bus::MIN_VOLTAGE_VIOLATION) {
			continue;
		}

		if (v == Bus::MIN_VOLTAGE_GR
				&& modBus->GetStatus() == Bus::MAX_VOLTAGE_VIOLATION) {
			continue;
		}
		std::cout << "Bus = " << crtBus->GetBus().m_nin << ", Mod = "
				<< modBus->GetStatus() << ", V = " << v << ", Vt = "
				<< vt(i) << std::endl;
		if (vt(i) > maxQ || maxBus == 0) {
			maxQ = vt(i);
			maxBus = i + 1;
		}
	}
	std::cout << "MaxBus = " << maxBus << std::endl;
	return maxBus;
}
}
