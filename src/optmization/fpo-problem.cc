#include "fpo-problem.h"
#include <vector>

using namespace arma;
using namespace std;
using namespace load;
using namespace model;

namespace opt {

FpoProblem::FpoProblem() {
	// TODO Auto-generated constructor stub
	m_lf = new LoadFlow();
	m_w = 0;
	m_nVar = 0;
	m_type = MIN;
}

FpoProblem::FpoProblem(string cdf, double w) {
	// TODO Auto-generated constructor stub
	m_lf = new LoadFlow();
	m_lf->Prepare(cdf);
	m_cdf = cdf;
	m_w = w;
	m_nVar = 0;
	m_type = MIN;
}

FpoProblem::~FpoProblem() {
	// TODO Auto-generated destructor stub
	m_lf = 0;
	m_cdf.clear();
}

FpoProblem::FpoProblem(LoadFlow* lf, double w) {
	m_lf = lf;
	m_w = w;
	m_nVar = 0;
	m_type = MIN;
}

double FpoProblem::CostFunction(Particle* p) {
	Graph* graph = m_lf->GetGraph();
	if(p != 0) {
		vec pos = p->GetPosition();
		for(int i = 0; i < graph->GetNumBus(); i++) {
			Bus* bus = graph->GetBus(i + 1);
			if(bus->GetType() == Bus::SLACK || bus->GetType() == Bus::GENERATION) {
				bus->SetVCalc(pos.at(bus->GetOrdG()));
			}
		}
		m_lf->Execute();
	}
	double fitness = 0;

	std::vector<Aresta_t> lista;
	int numBus = (int) graph->GetNumBus();

	for (int k = 0; k < numBus; k++) {
		Bus* busK = graph->GetBus(k + 1);
		//list.insert(k + 1);
		std::vector<Branch*> branches = busK->GetBranches();
		std::vector<Bus*> neighbors = busK->GetNeighbors();
		int numBranch = (int) branches.size();

		for (int m = 0; m < numBranch; m++) {
			Bus* busM = neighbors.at(m);
			Aresta_t list0(busK->GetBus().m_nin, busM->GetBus().m_nin);
			Aresta_t list1(busM->GetBus().m_nin, busK->GetBus().m_nin);
			bool flag = false;
			for (int i = 0; i < lista.size(); i++) {
				if (lista.at(i).equals(list0) || lista.at(i).equals(list1)) {
					flag = true;
					break;
				}
			}
			if (flag) {
				continue;
			}
			lista.push_back(list0);
			lista.push_back(list1);

			DBranch_t dataBranch = branches.at(m)->GetBranch();
			// TRANSMISSION_LINE = 0
			if (dataBranch.m_tipo == 0) {
				double theta_km = busK->GetACalc() - busM->GetACalc();
				fitness += dataBranch.m_g * (pow(busK->GetVCalc(), 2)
						   + pow(busM->GetVCalc(), 2) - 2 * busK->GetVCalc()
						   * busM->GetVCalc() * cos(theta_km));
			}
		}
	}
	// Penalização
	for (int i = 0; i < numBus; i++) {
		Bus* bus = graph->GetBus(i + 1);
		DBus_t dataBus = bus->GetBus();
		if (bus->GetType() == Bus::LOAD && bus->GetVCalc() < dataBus.m_vmin) {
			fitness = fitness + m_w * pow((bus->GetVCalc() - dataBus.m_vmin), 2);
		} else if (bus->GetType() == Bus::LOAD && bus->GetVCalc() > dataBus.m_vmax) {
			fitness = fitness + m_w * pow((bus->GetVCalc() - dataBus.m_vmax), 2);
		}
	}

	lista.clear();
	return fitness;
}

void FpoProblem::SetW(double w) {
	m_w = w;
}

void FpoProblem::SetCdf(string cdf) {
	m_cdf = cdf;
}

double FpoProblem::GetNVar() const {
	return m_nVar;
}

void FpoProblem::SetNVar(double nVar) {
	m_nVar = nVar;
}

vec FpoProblem::GetVarMax() const {
	return m_varMax;
}

void FpoProblem::SetVarMax(vec varMax) {
	m_varMax = varMax;
}

vec FpoProblem::GetVarMin() const {
	return m_varMin;
}

void FpoProblem::SetVarMin(vec varMin) {
	m_varMin = varMin;
}

vec FpoProblem::GetVarSize() const {
	return m_varSize;
}

void FpoProblem::SetVarSize(vec varSize) {
	m_varSize = varSize;
}

vec FpoProblem::GetMaxVelocity(void) const {
	return m_maxVelocity;
}

vec FpoProblem::GetMinVelocity(void) const {
	return m_minVelocity;
}

TypeOpt FpoProblem::GetType(void) const {
	return m_type;
}

void FpoProblem::SetType(TypeOpt type) {
	m_type = type;
}

VLim FpoProblem::CalcVelLimits(void) {
	m_maxVelocity = 0.2*(m_varMax - m_varMin);
	m_minVelocity = -m_maxVelocity;
	VLim lim(m_minVelocity, m_maxVelocity);
	return lim;
}

bool FpoProblem::Compare(double value1, double value2) {
	switch (m_type) {
		case MIN:
			return value1 < value2;
			break;
		case MAX:
			return value1 > value2;
			break;
		default:
			return value1 > value2;
			break;
	}
}

void FpoProblem::Print(void) {
	Graph* graph = m_lf->GetGraph();
	int numBus = (int) graph->GetNumBus();

	for (int k = 0; k < numBus; k++) {
		Bus* busK = graph->GetBus(k + 1);
		busK->Print();
	}
}

} /* namespace opt */
