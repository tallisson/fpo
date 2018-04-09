/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "fpo.h"
#include "report.h"
#include "../load/calc.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace model;
using namespace std;

namespace fpo {

Fpo::Fpo() {
	// TODO Auto-generated constructor stub
	m_erro = 0.001;
	m_maxIter = 100;
	m_w = 85;
	m_ro = 1;
	m_co = 0.00124;
	m_k = 0;
	m_verbose = false;
	lf = new LoadFlow();
	m_calc = new Calc();
	m_foMod = zeros(m_maxIter);
	m_perdas = zeros(m_maxIter);
}

Fpo::~Fpo() {
	m_erro = 0;
	m_maxIter = 0;
	m_w = 0;
	m_ro = 0;
	m_co = 0;

	m_u.clear();
	m_u.clear();
	m_x.clear();
	m_foMod.clear();
	m_perdas.clear();

	delete lf;
	lf = 0;

	delete m_calc;
	m_calc = 0;
}

void Fpo::SetErro(double erro) {
	m_erro = erro;
}

double Fpo::GetErro(void) const {
	return m_erro;
}

void Fpo::SetW(double w) {
	m_w = w;
}

double Fpo::GetW(void) const {
	return m_w;
}

void Fpo::SetMaxIter(double maxIter) {
	m_maxIter = maxIter;
}

double Fpo::GetMaxIter(void) const {
	return m_maxIter;
}

void Fpo::SetRo(double ro) {
	m_ro = ro;
}

double Fpo::GetRo(void) const {
	return m_ro;
}

void Fpo::SetVerbose(bool verbose) {
	m_verbose = verbose;
}

bool Fpo::GetVerbose(void) const {
	return m_verbose;
}

void Fpo::SetU(vec u) {
	m_u = u;
}

vec Fpo::GetU(void) const {
	return m_u;
}

void Fpo::Execute(std::string cdf) {
	LoadFlow* lf = new LoadFlow;
	lf->Prepare(cdf);
	Graph* graph = lf->GetGraph();
	Bus* slack = graph->GetSlack();

	if (m_verbose) {
		printf("w = %.4f \n", m_w);
	}

	m_x = zeros(graph->GetNumPV() + 1);

	if (m_verbose) {
		printf("-----------------------------------------\n");
		printf("Iteração: \n");

		printf("\nVariáveis de controle \n");
		printf("Barra V\n");
		printf("[k] [p.u.]\n");
		int size = (int) graph->GetBuses().size();
		for (int i = 1; i <= size; i++) {
			Bus* bus = graph->GetBus(i);
			if (bus->GetType() == Bus::GENERATION) {
				printf("%2d %12.4f\n", bus->GetBus().m_nin, bus->GetVCalc());
			}
		}
	}

	int size = (int) graph->GetBuses().size();
	for (int i = 1; i <= size; i++) {
		Bus* bus = graph->GetBus(i);
		if (bus->GetOrdG() != -1) {
			bus->AddControl(bus->GetVCalc());
			m_x(bus->GetOrdG()) = bus->GetVCalc();
		}
	}

	int flag = 0;
	int k = 0;
	while (flag == 0) {
		if (m_verbose) {
			printf("\n\nVariáveis dependentes: ");
			m_x.print(cout);
		}

		int conv = lf->Execute();
		if (m_verbose) {
			for (int i = 1; i <= size; i++) {
				Bus* bus = graph->GetBus(i);
				printf("%.10f\t", bus->GetVCalc());
			}

			printf("\n");
		}

		if (conv == 1) {
			lf->GetJac()->Zeros(graph);
			mat jac = lf->GetJac()->CalcJac(graph);
			slack->CalcPG();
			if (m_verbose) {
				cout << "Derivada da f.o. modificada em relação a: " << m_x
						<< endl;
			}
			vec dfdx = m_calc->GrafX(graph, m_w);
			vec lamb = inv(trans(jac)) * -dfdx;

			if (m_verbose) {
				cout << "Cálculo do gradiente reduzido:\n";
			}

			vec dLdu = m_calc->GradLU(graph, lamb);

			double fo = m_calc->Fitness(graph, 0);
			//cout << "Perdas = " << fo << endl;
			if (m_verbose) {
				cout << "Perdas = " << fo << endl;
			}
			m_foMod(k) = m_calc->Fitness(graph, m_w);
			m_perdas(k) = fo;
			//break;
			int nElem = (int) dLdu.n_rows;
			vec gradRed = zeros<vec>(nElem);
			gradRed = dLdu;

			int size = (int) graph->GetNumBus();
			for (int i = 0; i < size; i++) {
				Bus* bus = graph->GetBus(i + 1);
				if (bus->GetType() != Bus::LOAD
						&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
					if (bus->GetVCalc() - m_erro <= bus->GetBus().m_vmin) {
						if (gradRed(bus->GetOrdG()) >= 0) {
							gradRed(bus->GetOrdG()) = 0;
						}
					} else if (bus->GetVCalc() + m_erro
							>= bus->GetBus().m_vmax) {
						if (gradRed(bus->GetOrdG()) <= 0) {
							gradRed(bus->GetOrdG()) = 0;
						}
					}
				}
			}

			if (m_verbose) {
				cout << "Teste de convergência: " << endl;
			}
			ostringstream os;
			os.precision(13);
			os.setf(ios::fixed);
			os << norm(gradRed);
			double normV = atof(os.str().c_str());
			cout << "Norm = " << normV << endl;

			if (normV <= m_erro) {
				conv = 1;
			} else {
				conv = 0;
			}

			if (conv == 1) {
				flag = 1;
			} else if (k != m_maxIter) {
				k = k + 1;
				m_w = m_w * m_ro;
				//cout << "m_w = " << m_w << endl;

				double c = m_co;
				int numB = (int) graph->GetNumBus();
				Data_t bu = m_calc->Busca(lf, graph, m_w, c, m_erro, dLdu);
				cout << bu.print() << endl;
				if (bu.m_conv == 1) {
					for (int i = 0; i < numB; i++) {
						Bus* bus = graph->GetBus(i + 1);
						if (bus->GetType() != Bus::LOAD
								&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
							bus->SetVCalc(
									bus->GetVPreBusca()
											- bu.m_c * dLdu(bus->GetOrdG()));
						} else {
							bus->SetVCalc(bus->GetVPreBusca());
						}
						bus->SetACalc(bus->GetAPreBusca());
						bus->Print();
					}

					for (int i = 0; i < numB; i++) {
						Bus* bus = graph->GetBus(i + 1);
						if (bus->GetType() != Bus::LOAD
								&& bus->GetType() != Bus::LOSS_CONTROL_REACT
								&& bus->GetVCalc() < bus->GetBus().m_vmin) {
							bus->SetVCalc(bus->GetBus().m_vmin);
							if (m_verbose) {
								printf("\nV%f < V_min! -> V%f = V_min",
										bus->GetVCalc(), bus->GetVCalc());
							}
						} else if (bus->GetType() != Bus::LOAD
								&& bus->GetType() != Bus::LOSS_CONTROL_REACT
								&& bus->GetVCalc() > bus->GetBus().m_vmax) {
							bus->SetVCalc(bus->GetBus().m_vmax);
							if (m_verbose) {
								printf("\nV%f > V_max! -> V%f = V_max",
										bus->GetVCalc(), bus->GetVCalc());
							}
						}
					}
					if (m_verbose) {
						printf("\n\nVariáveis de controle 'u':");
						printf("\nBarra    V");
						printf("[k]      [p.u.]");
						for (int i = 0; size; i++) {
							Bus* bus = graph->GetBus(i + 1);
							if (bus->GetType() != Bus::GENERATION) {
								printf("%2d %12.4f", bus->GetBus().m_nin,
										bus->GetVCalc());
							}
						}
					}

					for (int i = 0; i < numB; i++) {
						Bus* bus = graph->GetBus(i + 1);
						if (bus->GetType() != Bus::GENERATION) {
							bus->AddControl(bus->GetVCalc());
						}
					}

					if (m_verbose) {
						printf("\n\n-----------------------------------------");
						printf("\nIteração %d:", k);
					}

				} else {
					flag = 4;
				}
			} else {
				flag = 2;
			}

		} else {
			flag = 3;
		}
		printf("Flag = %d AND k = %d", flag, k);
	}

// Relatório de saída:
	m_verbose = true;
	if (m_verbose) {
		if (flag == 1) {
			printf(
					"\n\nO método do Gradiente Reduzido convergiu em %d iterações com E = %.4f",
					k, m_erro);
		} else {
			if (flag == 2) {
				printf(
						"\n\nO número máximo de iterações do método do Gradiente Reduzido foi atingido...");
			} else if (flag == 3) {
				printf(
						"\n\nO processo foi interrompido porque o problema de Fluxo de Carga (Passo 2) não convergiu...");
			} else {
				printf(
						"\n\nO processo foi interrompido porque a busca unidimensional não conseguiu determinar um passo descendente...");
			}
		}

		if (flag != 3) {
			Report(lf, k);
		}
	}
}

void Fpo::Report(LoadFlow* lf, int numIt) {
	printf("\n\nResumo do Relatório de Convergência:");
	printf("\n\nRelatório de Saída:");
	printf("\nIteração   V1        V2        FO Mod    Perdas\n");
	printf("[it]       [p.u.]    [p.u.]    [p.u.]    [p.u.]\n");
	for (int i = 0; i <= numIt; i++) {
	  printf("%3d %13.4f %9.4f %9.4f %9.4f\n", i, 0.0, 0.0, m_foMod.at(i),m_perdas.at(i));
	}
	printf("\n\nRelatório de Saída:\n");
	int numB = (int) lf->GetGraph()->GetNumBus();
	for (int i = 0; i < numB; i++) {
		Bus* bus = lf->GetGraph()->GetBus(i + 1);
		bus->Print();
	}

	// Geração de potência ativa:
	Bus* slack = lf->GetGraph()->GetSlack();
	slack->CalcPG();

	// Geração de potência reativa:
	for (int i = 0; i < numB; i++) {
		Bus* bus = lf->GetGraph()->GetBus(i + 1);

		if (bus->GetType() != Bus::LOAD
				|| bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			bus->CalcQg();
		}
	}

	// Impressão na tela do relatório de saída:

	printf(
			"\nBarra   V         ang       ang       Pg        Qg        Pc        Qc\n");
	printf(
			"[k]     [p.u.]    [graus]   [rad]     [p.u.]    [p.u.]    [p.u.]    [p.u.]\n");
	for (int i = 0; i < numB; i++) {
		Bus* bus = lf->GetGraph()->GetBus(i + 1);
		printf("%2d %11.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
				bus->GetBus().m_nin, bus->GetVCalc(),
				bus->GetACalc() * (180 / M_PI), bus->GetACalc(),
				bus->GetPCalc(), bus->GetQCalc(), bus->GetBus().m_pc,
				bus->GetBus().m_qc);
	}

	// Cálculo dos fluxos nas linhas de transmissão e perdas ativas:
	DataLoss_t data = lf->CalcLosses();
	vec pkm = data.m_pkm;
	int numD = (int) pkm.size();
	printf("\n\nPkm       Pmk       Perdas    Qkm       Qmk");
	printf("\n[p.u.]    [p.u.]    [p.u.]    [p.u.]    [p.u.]\n");
	for (int i = 0; i < numD; i++) {
		printf("%6.4f %9.4f %9.4f %9.4f %9.4f\n", pkm.at(i), data.m_pmk.at(i),
				data.m_l.at(i), data.m_qkm.at(i), data.m_qmk.at(i));
	}

	printf("Perdas Totais = %.4f [MW]", data.m_total);
}

}
