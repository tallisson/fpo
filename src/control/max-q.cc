/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "max-q.h"

#include <armadillo>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

using namespace arma;

namespace control {

MaxQ::MaxQ() {
}

MaxQ::~MaxQ() {
}

bool MaxQ::DoControl(mat jqv, Graph* graph) {
	bool control = false;
	vec deltaV = zeros<vec>(jqv.n_cols);
	for (int i = 1; i <= graph->GetNumBus(); i++) {
		Bus* bus = graph->GetBus(i);
		if (bus->GetType() != Bus::LOAD
				&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			continue;
		}
		double dsv = bus->CalcDsv();
		if (dsv != 0) {
			deltaV(i - 1) = dsv;
			control = true;
		}
	}
	std::cout.precision(11);
	std::cout.setf(ios::fixed);

	deltaV.raw_print(cout, "Delta V");
	std::cout << std::endl;
	if (control == true) {
		Bus* maxVlt = MaxDsv(graph);
		vec deltaQ = (jqv * deltaV);
		deltaQ.raw_print(cout, "Delta Q:");
		std::cout << std::endl;
		uint32_t idBus = MaxV(graph, deltaQ, maxVlt);
		if (idBus == 0) {
			std::cout << "Error definition bus MaxQ" << std::endl;
			return false;
		}

		std::cout << "Max Q = " << deltaQ(idBus - 1) << std::endl;
		vec auxQ = zeros<vec>(deltaQ.n_elem);
		auxQ(idBus - 1) = deltaQ(idBus - 1);

		jqv.raw_print(cout, "jqv:");
		std::cout << std::endl;
		mat iv = inv(jqv);
		iv.raw_print(cout, "jqv^-1:");
		std::cout << std::endl;

		vec deltaVIjt = inv(jqv) * auxQ;

		double m_alpha = 1;
		double value = fabs(deltaVIjt(idBus - 1));
		while (m_alpha * value < LIMIAR) {
			m_alpha++;
		}
		//std::cout << "Variação de Tensão => " << value << ", alpha = " << m_alpha << "\n";

		Bus* bus = graph->GetBus(idBus);
		value = deltaVIjt(idBus - 1);
		if (maxVlt->GetStatus() == Bus::MIN_VOLTAGE_VIOLATION && value < 0) {
			value = fabs(value);
		}

		if (maxVlt->GetStatus() == Bus::MAX_VOLTAGE_VIOLATION && value > 0) {
			value *= -1;
		}
		double v;
		v = bus->GetVCalc();
		double newValue = v + (value * m_alpha);
		if (newValue < Bus::MIN_VOLTAGE_GR) {
			newValue = Bus::MIN_VOLTAGE_GR;
		}
		if (newValue > Bus::MAX_VOLTAGE_ONS) {
			newValue = Bus::MAX_VOLTAGE_ONS;
		}

		std::cout << "Bus => " << bus->GetBus().m_nin << ", Voltage + value = "
				<< (newValue) << std::endl;
		bus->SetVCalc(newValue);
	}
	return control;
}

}
