/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "vsf.h"
#include "../model/graph.h"

#include <armadillo>
#include <map>
#include <math.h>
#include <vector>

using namespace arma;

namespace control {

Vsf::Vsf() {
}

Vsf::~Vsf() {
}

bool Vsf::DoControl(mat jqv, Graph* graph) {
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
		mat invJqv = inv(jqv);

		jqv.raw_print(cout, "jqv:");
		std::cout << std::endl;
		invJqv.raw_print(cout, "jqv^-1:");

		subview_col<double> s = invJqv.col(maxVlt->GetBus().m_nin - 1);
		vec vsf = zeros<vec>(s.n_elem);
		for (uint32_t i = 0; i < s.n_elem; i++) {
			vsf(i) = s(i);
		}

		uint32_t idBus = MaxV(graph, vsf, maxVlt);
		if (idBus == 0) {
			std::cout << "Error definition bus Vsf" << std::endl;
			return false;
		}
		double t;
		t = maxVlt->GetVCalc();
		std::cout << "Max Vlt = " << maxVlt->GetBus().m_nin << " = " << t
				<< ", Max = " << vsf(idBus - 1) << ", IdBus = " << idBus
				<< std::endl;
		vec aux = zeros<vec>(vsf.n_elem);
		aux(idBus - 1) = vsf(idBus - 1);
		vec deltaVIjt = inv(jqv) * aux;

		double m_alpha = 1;
		double value = fabs(deltaVIjt(idBus - 1));
		while (m_alpha * value < LIMIAR) {
			m_alpha++;
		}
		std::cout << "Variação de Tensão => " << value << ", alpha = "
				<< m_alpha << "\n";

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

		std::cout << "Value = " << v << " Bus = " << bus->GetType()
				<< ", Voltage + value = " << (newValue) << std::endl;
		bus->SetVCalc(newValue);
	}

	return control;
}

}
