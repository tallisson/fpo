/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "load-flow.h"

#include "../control/q-control.h"

#include <iostream>
#include <armadillo>
#include <map>

#include <sstream>

using namespace arma;

namespace load {

LoadFlow::LoadFlow(void) {
	m_precision = 0.00001;
	m_iter = 0;
	m_maxIter = 10;
	m_qControl = new QControl();
	m_graph = new Graph();
	m_mismatches = new Mismatch();
	m_jac = new Jacobian();
	m_verbose = false;
	m_hasQControl = false;
	m_totalL = 0.0;
	m_vControl = NULL;
}

LoadFlow::~LoadFlow(void) {
	delete m_graph;
	m_graph = NULL;
	delete m_qControl;
	m_qControl = NULL;
	delete m_mismatches;
	m_mismatches = NULL;
	delete m_jac;
	m_jac = NULL;
}

Graph* LoadFlow::GetGraph(void) const {
	return m_graph;
}

void LoadFlow::SetQControl(bool hasQControl) {
	m_hasQControl = hasQControl;
}

Control* LoadFlow::GetQControl(void) const {
	return m_qControl;
}

void LoadFlow::SetQControl(Control* qControl) {
	m_qControl = qControl;
}

void LoadFlow::InitX0(void) {
	// Passo 2 - Arbitrar um ponto inicial:
	for (int i = 0; i < m_graph->GetNumBus(); i++) {
		Bus* bus = m_graph->GetBus(i + 1);
		if (bus->GetType() != Bus::SLACK) {
			/*double angSlack =
					m_graph->GetBus(m_graph->GetPosSlack())->GetBus().m_ang;*/
			double angSlack = 0;
			bus->SetACalc(angSlack);
		}

		if (bus->GetType() == Bus::LOAD) {
			bus->SetVCalc(1.0);
		}
	}
}

void LoadFlow::Prepare(std::string cdf) {
	// Definir topologia completa da rede.
	Utils* utils = new Utils();
	m_sts = utils->Read(cdf);

	m_graph->SetPosSlack(m_sts.m_posSlack);
	if(m_verbose) {
		cout << "Dados das Barras " << endl;
	}
	int size = (int)m_sts.buses.size();
	for (int i = 0; i < size; i++) {
		DBus_t busData = m_sts.buses.at(i);
		Bus* bus = new Bus();
		bus->SetBus(busData);
		m_graph->AddBus(bus);
		if(m_verbose) {
			bus->Print();
		}
	}
	if(m_verbose) {
		cout << "Dados dos Branches " << size << std::endl;
	}
	size = m_sts.branches.size();
	for (int i = 0; i < size; i++) {
		DBranch_t dataBranch = m_sts.branches.at(i);
		Branch* branch = new Branch();
		branch->SetBranch(dataBranch);
		m_graph->Assoc(branch);
		if (dataBranch.m_tipo == 0) {
			dataBranch.m_tap = 1.0;
			branch->SetBranch(dataBranch);

			Bus* busNi = m_graph->GetBus(dataBranch.m_ni);
			Bus* busNf = m_graph->GetBus(dataBranch.m_nf);
			busNi->SetTap(Bus::IMP);
			busNf->SetTap(Bus::IMP);
		} else if (dataBranch.m_tipo == 1) {
			Bus* busNi = m_graph->GetBus(dataBranch.m_ni);
			Bus* busNf = m_graph->GetBus(dataBranch.m_nf);
			busNi->SetTap(Bus::TAP);
			busNf->SetTap(Bus::IMP);
		}
		if(m_verbose) {
			branch->Print();
		}
	}
}

void LoadFlow::InitJ(void) {
	int nPQ = m_graph->GetNumPQ();
	int nPV = m_graph->GetNumPV();
	int size = m_graph->GetNumPQ() * 2 + m_graph->GetNumPV();
	int numB = m_graph->GetNumBus();
	m_jac->SetMatrix(size, size);

	m_jac->SetJ1((nPQ + nPV), (nPQ + nPV));
	m_jac->SetJ2((nPQ + nPV), numB);
	m_jac->SetJ3(numB, (nPQ + nPV));
	m_jac->SetJ4(numB, numB);

}

void LoadFlow::Reset(void) {
	int size = m_graph->GetNumPQ() * 2 + m_graph->GetNumPV();
	m_jac->Zeros(size, size);
	m_b = zeros<vec>(size);
}

int LoadFlow::Execute() {
	/*
	 * Calcular o vetor dos mismatches, a matriz Jacobiana e resolver
	 * o sistema de equações lineares:
	 */
	m_conv = 1;
	InitX0();
	InitJ();

	m_b = m_mismatches->CalcMismatches(m_graph);
	std::cout.precision(11);
	std::cout.setf(ios::fixed);
	if(m_verbose) {
		m_b.raw_print(cout, "\n160-lf.cc Erro: ");
		cout << "\n" << endl;
	}
	bool execute = false;
	int nextIter, nextCrt;
	nextCrt = 0;
	if(m_verbose) {
		std::cout << "Estado inicial" << std::endl;

		for (int i = 1; i <= m_graph->GetNumBus(); i++) {
			Bus* bus = m_graph->GetBus(i);
			bus->Print();
		}
		m_b.raw_print(cout, "Erros:");
		std::cout << std::endl;
		std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
	do //while (execute)
	{
		nextIter = 0;
		m_iter = 0;
		while (nextIter == 0) {
			mat m = m_jac->CalcJac(m_graph);
			vec dx = m_jac->SolveSys(m_b);
			if(m_verbose) {
				m.raw_print(cout, "Jacobiana: ");
				std::cout << std::endl;
				std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
				m_b.raw_print(cout, "Erros:");
				std::cout << std::endl;
				std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
				dx.raw_print(cout, "X:");
				std::cout << std::endl;
				std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
			}

			// Atualizar 'a' e 'V':
			for (int i = 0; i < m_graph->GetNumBus(); i++) {
				Bus* bus = m_graph->GetBus(i + 1);
				if (bus->GetType() != Bus::SLACK) {
					double angD = bus->GetACalc();

					angD += dx(bus->GetBus().m_ord);
					bus->SetACalc(angD);
				}

				if (bus->GetType() != Bus::SLACK
						&& bus->GetType() != Bus::GENERATION) {
					double vD = bus->GetVCalc();

					int ind = m_graph->GetNumBus() - 1
							+ bus->GetBus().m_ordPQ;
					vD += dx(ind);
					bus->SetVCalc(vD);
				}
			}
			bool crt = false;
			if (m_hasQControl == true) {
				crt = m_qControl->DoRestore(m_graph);
				if (crt == true) {
					InitJ();
				}
			}

			m_iter++;

			// Qlim
			if (m_hasQControl == true) {
				crt = m_qControl->DoControl(m_graph);
				if (crt == true) {
					InitJ();
				}
			}

			/*
			 * Cálculo do vetor dos mismatches com os valores corrigidos de
			 * 'a' e 'V':
			 */
			m_b = m_mismatches->CalcMismatches(m_graph);
			//m_b.raw_print(cout, "(238 lf.cc) Erro: ");

			if(m_verbose) {
				std::cout << "Erros: \n" << m_b;
				std::cout << "+++++++++++++++++++++++++++++++++++++++++\n";

				for (int i = 0; i < m_graph->GetNumBus(); i++) {
					Bus* bus = m_graph->GetBus(i + 1);
					bus->Print();
				}
			}
			// Teste de convergência:
			double maxB = max(abs(m_b));
			if (maxB <= m_precision) {
				nextIter = 1;
				m_conv = 1;
			} else {
				m_jac->Zeros();
				int nPQ = m_graph->GetNumPQ();
				int nPV = m_graph->GetNumPV();
				int numB = m_graph->GetNumBus();
				m_jac->SetJ1((nPQ + nPV), (nPQ + nPV));
				m_jac->SetJ2((nPQ + nPV), numB);
				m_jac->SetJ3(numB, (nPQ + nPV));
				m_jac->SetJ4(numB, numB);
				nextIter = 0;

				// Critério de saída do laço:
				if (m_iter == m_maxIter) {
					nextIter = 2;
				}
				m_conv = 0;
			}
			if (nextIter != 0) {
				CalcLosses();
			}
		}
		if (m_vControl != NULL) {
			execute = false;
			if (nextCrt < 2) {
				execute = m_vControl->DoControl(m_jac->GetJqv(), m_graph);
			}
		}
		if (execute == true) {
			nextCrt++;
		}
	} while (execute);
	if (nextIter == 1 && m_verbose) {
		std::cout << "O método de Newton-Raphson convergiu em " << m_iter
				<< " iterações" << std::endl;
	} else if(m_verbose){
		std::cout
				<< "O número máximo de iterações foi atingido e o método de Newton-Raphson não convergiu..."
				<< std::endl;
	}
	m_verbose = true;
	if(m_verbose) {
		for (int i = 0; i < m_graph->GetNumBus(); i++) {
			Bus* bus = m_graph->GetBus(i + 1);
			bus->Print();
		}
	}

	return m_conv;
}

void LoadFlow::CalcLosses(void) {
	std::vector<double> p_km, p_mk, q_km, q_mk, pL, losses;
	m_totalL = 0;
	std::vector<Branch* > branches = m_graph->GetBranches();
	for (int i = 0; i < m_graph->GetNumBranch(); i++) {
		Branch* branch = branches.at(i);
		Bus* busK = m_graph->GetBus(branch->GetBranch().m_ni);
		Bus* busM = m_graph->GetBus(branch->GetBranch().m_nf);

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
		m_totalL += l;
	}

	m_totalL *= m_sts.m_baseMVA;
	if (m_verbose) {
		//@ToDo
	}
	if(m_verbose) {
		std::cout << "Perda Total: " << m_totalL << std::endl;
	}
}

void LoadFlow::SetPrecision(double precision) {
	m_precision = precision;
}

void LoadFlow::SetDir(std::string dir) {
	m_dir = dir;
}

void LoadFlow::SetFile(std::string file) {
	m_file = file;
}

void LoadFlow::SetVControl(VControl* vControl) {
	m_vControl = vControl;
}

void LoadFlow::SetVerbose(bool verbose) {
	m_verbose = verbose;
}

bool LoadFlow::GetVerbose(void) const {
	return m_verbose;
}

mat LoadFlow::ExecuteJ(void) {
	return m_jac->CalcJac(m_graph);
}

int LoadFlow::GetConv(void) const {
	return m_conv;
}

Jacobian* LoadFlow::GetJac(void) const {
	return m_jac;
}

}

