/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "calc.h"
#include "jacobian-u.h"

#include <iostream>
#include <stdio.h>
#include <set>

namespace load {

Calc::Calc() {
	// TODO Auto-generated constructor stub
	m_verbose = false;
}

Calc::~Calc() {
	// TODO Auto-generated destructor stub
}

void Calc::IsVerbose(bool verbose) {
	m_verbose = verbose;
}

// Determinação do vetor gradiente da função objetivo em relação às variáveis dependentes 'x':
vec Calc::GrafX(Graph* graph, double w) {
	vec y = zeros(graph->GetNumPV() + 2 * graph->GetNumPQ());
	int size = (int) graph->GetBuses().size();
	//  Derivada da função objetivo em relação a 'ak':
	for (int k = 0; k < size; k++) {
		Bus* busK = graph->GetBus(k + 1);
		if (busK->GetType() != Bus::SLACK) {
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int numBranch = (int) branches.size();
			for (int m = 0; m < numBranch; m++) {
				Bus* busM = neighbors.at(m);
				Branch* branch = branches.at(m);
				DBranch_t dataBranch = branch->GetBranch();
				double theta_km = busK->GetACalc() - busM->GetACalc();
				// função objetivo em relação a 'ak':
				// y(k-1) = y(k-1) + 2*g(km)*V(k)*V(m)*sin(akm);
				y(busK->GetBus().m_ord) += 2 * dataBranch.m_g * busK->GetVCalc()
						* busM->GetVCalc() * sin(theta_km);
			}
		}
	}

	// Derivada da função objetivo modificada em relação a 'Vk':
	for (int k = 0; k < size; k++) {
		Bus* busK = graph->GetBus(k + 1);
		if (busK->GetType() == Bus::LOAD) {
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int numBranch = (int) branches.size();
			for (int m = 0; m < numBranch; m++) {
				Bus* busM = neighbors.at(m);
				Branch* branch = branches.at(m);
				DBranch_t dataBranch = branch->GetBranch();
				double theta_km = busK->GetACalc() - busM->GetACalc();
				// função objetivo em relaçao a 'Vk':
				// y(nb-1+ordPQ(k)) = y(nb-1+ordPQ(k)) + 2*g(km)*(V(k) - V(m)*cos(akm));
				y(graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ) += 2
						* dataBranch.m_g
						* (busK->GetVCalc() - busM->GetVCalc() * cos(theta_km));
			}

			// derivada da função penalidade em relação a 'Vk':
			if (busK->GetVCalc() < busK->GetBus().m_vmin) {
				y(graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ) += 2 * w
						* (busK->GetVCalc() - busK->GetBus().m_vmin);
				//y(nb-1+ordPQ(k)) = y(nb-1+ordPQ(k)) + 2*w*(V(k)-V_min(k));
			} else if (busK->GetVCalc() > busK->GetBus().m_vmax) {
				y(graph->GetNumBus() - 1 + busK->GetBus().m_ordPQ) += 2 * w
						* (busK->GetVCalc() - busK->GetBus().m_vmax);
				//y(nb-1+ordPQ(k)) = y(nb-1+ordPQ(k)) + 2*w*(V(k)-V_max(k));
			}
		}
	}

	return y;
}

vec Calc::GradLU(Graph* graph, vec lamb) {
	// Determinação do vetor gradiente da função Lagrangiana em relação às variáveis 'u':

	// Cálculo do vetor gradiente da função objetivo em relação às variáveis 'u':
	vec dfdu = GradFU(graph);
	cout << "dfdu = \n" << dfdu << endl;

	// Cálculo da Jacobiana das equações do problema de FC em relação às variáveis 'u':
	JacobianU* jac = new JacobianU();
	mat jac_u = jac->CalcJac(graph);
	cout << "JacU \n" << jac_u << endl;

	vec y = dfdu + trans(jac_u) * lamb;
	//delete jac;
	return y;
}

vec Calc::GradFU(Graph* graph) {
	// Determinação do vetor gradiente da função objetivo em relação às variáveis de controle 'u':
	vec y = zeros(graph->GetNumPV() + 1);

	// Derivada da função objetivo em relação a 'Vk':
	int size = (int) graph->GetNumBus();
	for (int k = 0; k < size; k++) {
		Bus* busK = graph->GetBus(k + 1);
		if (busK->GetType() != Bus::LOAD) {
			std::vector<Branch*> branches = busK->GetBranches();
			std::vector<Bus*> neighbors = busK->GetNeighbors();
			int numBranch = (int) branches.size();
			for (int m = 0; m < numBranch; m++) {
				Bus* busM = neighbors.at(m);
				Branch* branch = branches.at(m);
				DBranch_t dataBranch = branch->GetBranch();
				double theta_km = busK->GetACalc() - busM->GetACalc();
				// função objetivo em relaçao a 'Vk':
				// y(ordG(k)) = y(ordG(k)) + 2*g(km)*(V(k) - V(m)*cos(akm));
				y(busK->GetOrdG()) += 2 * dataBranch.m_g
						* (busK->GetVCalc() - busM->GetVCalc() * cos(theta_km));
			}
		}
	}
	return y;
}

double Calc::Fitness(Graph* graph, double w) {
	// Cálculo da função objetivo:
	double fitness = 0;
	std::set<int> list;

	int numBus = (int) graph->GetNumBus();
	double v2 = graph->GetBus(2)->GetVCalc();
	for (int k = 0; k < numBus; k++) {
		Bus* busK = graph->GetBus(k + 1);
		list.insert(k + 1);
		std::vector<Branch*> branches = busK->GetBranches();
		std::vector<Bus*> neighbors = busK->GetNeighbors();
		int numBranch = (int) branches.size();

		for (int m = 0; m < numBranch; m++) {
			Bus* busM = neighbors.at(m);
			std::set<int>::iterator it = list.find(busM->GetBus().m_nin);
			if (it != list.end()) {
				continue;
			}
			DBranch_t dataBranch = branches.at(m)->GetBranch();

			// TRANSMISSION_LINE = 0
			if (dataBranch.m_tipo == 0) {
				double theta_km = busK->GetACalc() - busM->GetACalc();
				fitness = fitness + dataBranch.m_g * ( pow(busK->GetVCalc(), 2)
						  + pow(busM->GetVCalc(), 2) - 2 * busK->GetVCalc()* busM->GetVCalc() * cos(theta_km) );
				if(v2 > 1.10) {
					cout << "Fitness " << fitness << " AND " << busK->GetACalc() << endl;
				}
			}
		}
	}
	for (int i = 0; i < numBus; i++) {
		Bus* bus = graph->GetBus(i + 1);
		DBus_t dataBus = bus->GetBus();
		if (bus->GetType() == Bus::LOAD && bus->GetVCalc() < dataBus.m_vmin) {
			fitness = fitness + w * pow((bus->GetVCalc() - dataBus.m_vmin), 2);
		} else if (bus->GetType() == Bus::LOAD && bus->GetVCalc() > dataBus.m_vmax) {
			fitness = fitness + w * pow((bus->GetVCalc() - dataBus.m_vmax), 2);
		}
	}
	list.clear();

	return fitness;
}

Data_t Calc::Busca(LoadFlow* lf, Graph* graph, double w, double c, double erro,
		vec dLdu) {
	Data_t y;
	// Determinação do valor do parâmetro "c":
	double ro = 2;
	double min_c = 1e-6;
	double max_c = 5;
	double cont = 0;
	//Bus* slack = graph->GetSlack();

	// Cálculo da função objetivo modificada no ponto atual:
	double fo1 = this->Fitness(graph, w);

	if (m_verbose) {
		cout << "Busca Unidimensional: " << endl;
		cout << "fo1        c        fo2" << endl;
	}

	// Pré-busca
	int numBus = (int) graph->GetNumBus();
	for(int k = 0; k < numBus; k++) {
		Bus* bus = graph->GetBus(k + 1);
		bus->SetVPreBusca(bus->GetVCalc());
		bus->SetAPreBusca(bus->GetACalc());
	}

	vec u_o = zeros<vec>(dLdu.n_elem);
	for (int k = 0; k < numBus; k++) {
		Bus* bus = graph->GetBus(k + 1);
		if (bus->GetType() != Bus::LOAD) {
			u_o(bus->GetOrdG()) = bus->GetVCalc();
		}
	}
	//printf("c = %.12f, w = %.12f e erro = %.12f\n", c, w, erro);
	cout << "Busca = " << endl;

	// Determinação do valor do parâmetro "c":
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;
	//int cont = 0;

	while (flag1 == 0 || flag2 == 0) {
		// Ajuste das magnitudes de tensão das barras de geração:
		// (PV + slack)
		for (int k = 0; k < numBus; k++) {
			Bus* bus = graph->GetBus(k + 1);
			if (bus->GetType() != Bus::LOAD) {
				bus->SetVCalc(bus->GetVPreBusca() - c * dLdu(bus->GetOrdG()));
			}
		}

		// Correção dos limites violados:
		for (int k = 0; k < numBus; k++) {
			Bus* bus = graph->GetBus(k + 1);
			if (bus->GetType() != Bus::LOAD) {
				if (bus->GetVCalc() < bus->GetBus().m_vmin) {
					bus->SetVCalc(bus->GetBus().m_vmin);
				} else if (bus->GetVCalc() > bus->GetBus().m_vmax) {
					bus->SetVCalc(bus->GetBus().m_vmax);
				}
			}
			if(cont == 5) {
				printf("%.12f\n", bus->GetVCalc());
			}
		}

		// Cálculo do Fluxo de Carga ("atualização" das variáveis dependentes 'x'):
		lf->SetInit(true);
		int conv = lf->Execute();
		if (conv == 1) {
			double pgSlack = graph->GetSlack()->CalcPG();
			cout << "Slack Pg = " << pgSlack << endl;

			// Cálculo da função objetivo modificada no ponto candidato a ponto ótimo:
			double fo2 = Fitness(graph, w);
			cout << "Fo1 = " << fo1 << endl;
			cout << "Fo2 = " << fo2 << endl;
			cout << "Cont = " << cont << endl;
			m_verbose = true;
			if (m_verbose) {
				if (fo2 < fo1) {
					printf("\n%.6f   %6.4f   %.6f  <-  ok\n", fo1, c, fo2);
				} else {
					printf("\n%.6f   %6.4f   %.6f  <-   x\n", fo1, c, fo2);
				}
			}
			m_verbose = false;

			// Teste do passo descendente:
			if ((c > min_c) || (c < max_c)) {
				if (fo2 < fo1) {
					c *= ro;
					flag1 = 1;
				} else {
					c /= ro;
					flag2 = flag2 + 1;
				}
			} else {
				flag2 = 1;
				flag3 = 1;
			}
		} else {
			c /= ro;
			flag2 += 1;
		}

		if (m_verbose) {
			printf("\nc = %.4f\n", c);
		}
		cont++;
	}
	printf("C = %.12f, Flag 1 = %d, Flag 2 = %d, Flag 3 = %d\n", c, flag1, flag2, flag3);
	y.m_c = c;
	//return y;
	if (flag3 != 1 || flag1 == 1) {
		y.m_conv = 1;
	} else {
		y.m_conv = 0;
	}

	//printf("Conv = %.12f AND c = %.12f\n", y.m_conv, y.m_c);
	return y;
}

} /* namespace load */
