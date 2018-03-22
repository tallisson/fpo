/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "report.h"

#include <vector>

namespace fpo {

Report::Report() {
	// TODO Auto-generated constructor stub

}

Report::~Report() {
	// TODO Auto-generated destructor stub
}

void Report::Output(Graph* graph) {
	// Cálculo das grandezas desconhecidas associadas às barras slack e PV:

	// Cálculo dos fluxos nas linhas de transmissão e perdas ativas:
	int numBus = (int)graph->GetNumBus();
	printf("\nDados das barras\n");
	for(int i = 0; i < numBus; i++) {
		Bus* bus = graph->GetBus(i + 1);
		if(bus->GetType() == Bus::SLACK) {
			bus->CalcPg();
		}
		if(bus->GetType() != Bus::LOAD) {
			bus->CalcQg();
		}
		bus->Print();
	}
	std::vector<double> p_km, p_mk, q_km, q_mk, pL, losses;
	double totalL = 0;
	std::vector<Branch* > branches = graph->GetBranches();
	for (int i = 0; i < graph->GetNumBranch(); i++) {
		Branch* branch = branches.at(i);
		Bus* busK = graph->GetBus(branch->GetBranch().m_ni);
		Bus* busM = graph->GetBus(branch->GetBranch().m_nf);

		double vK, vM, aK, aM;
		vK = busK->GetVCalc();
		aK = busK->GetACalc();
		vM = busM->GetVCalc();
		aM = busM->GetACalc();

		p_km.push_back(branch->CalcPkmL(vK, vM, aK, aM));
		p_mk.push_back(branch->CalcPmkL(vK, vM, aK, aM));
		q_km.push_back(branch->CalcQkmL(vK, vM, aK, aM));
		q_mk.push_back(branch->CalcQmkL(vK, vM, aK, aM));

		double l = branch->CalcL(vK, vM, aK, aM);
		losses.push_back(l);
		totalL += l;
	}

	std::cout << "Perda Total: " << totalL << std::endl;
}

} /* namespace load */
