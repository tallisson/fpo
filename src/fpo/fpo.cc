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
	m_co = 0.000124;
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

/* Resolução do problema de Despacho Econômico com Representação da Rede
 * pelo método do Gradiente Reduzido (Dommel e Tinney):
 */
void Fpo::Execute(std::string cdf) {
	LoadFlow* lf = new LoadFlow;
	lf->Prepare(cdf);
	Graph* graph = lf->GetGraph();
	Bus* slack = graph->GetSlack();

	if (m_verbose) {
		printf("w = %.4f \n", m_w);
	}

	m_x = zeros(graph->GetNumPV() + 1);
	m_verbose = true;
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

	// Armazenamento das variáveis de controle a cada iteração: ctrls = zeros(max_iter, npv + 1);
	int size = (int) graph->GetBuses().size();
	for (int i = 1; i <= size; i++) {
		Bus* bus = graph->GetBus(i);
		if (bus->GetOrdG() != -1) {
			bus->AddControl(bus->GetVCalc());
			m_x(bus->GetOrdG()) = bus->GetVCalc();
		}
	}
	m_verbose = false;

	/* Declaração dos vetores para o armazenamento dos valores da função objetivo
	 modificada e da função perdas a cada iteração: fo_mod = zeros(max_iter, 1);
	 perdas = zeros(max_iter, 1); Início do processo iterativo: */
	int flag = 0;
	int k = 0;
	while (flag == 0) {
		// Passo 2 - Calcular as variáveis dependentes x em função de u:
		// (Resolução do problema de Fluxo de Carga)
		if (m_verbose) {
			printf("\n\nVariáveis dependentes: ");
			m_x.print(cout);
		}

		int conv = lf->Execute();
		m_verbose = true;
		if (m_verbose) {
			for (int i = 1; i <= size; i++) {
				Bus* bus = graph->GetBus(i);
				printf("%.10f\t", bus->GetVCalc());
			}
		}
		printf("\n");
		m_verbose = false;

		if (conv == 1) {
			/* Atualização da matriz Jacobiana para os novos valores de "V" e "a":
			 (Este passo é desnecessário quando o ponto de operação já está próximo ao ótimo.)*/
			//mat jac = lf->ExecuteJ();
			lf->GetJac()->Zeros(graph);
			mat jac = lf->GetJac()->CalcJac(graph);

			cout.precision(12);
			cout.setf(ios::fixed);
			jac.raw_print(cout, "Jacobiana ");
			slack->CalcPG();
			if (m_verbose) {
				cout << "Derivada da f.o. modificada em relação a: " << m_x
						<< endl;
			}
			vec dfdx = m_calc->GrafX(graph, m_w);
			/*printf("dfdx: \n");
			 for(int t = 0; t < dfdx.size(); t++) {
			 printf("%.10f\n", dfdx(t));
			 }*/
			vec lamb = inv(trans(jac)) * -dfdx;
			/*printf("Lamb: \n");
			 for(int t = 0; t < lamb.size(); t++) {
			 printf("%.10f\n", lamb(t));
			 }*/

			if (m_verbose) {
				cout << "Cálculo do gradiente reduzido:\n";
			}

			vec dLdu = m_calc->GradLU(graph, lamb);
			dLdu.raw_print(cout, "dLdu ");

			// Cálculo da função objetivo no ponto atual (sem a função penalidade -> w=0):
			double fo = m_calc->Fitness(graph, 0);
			cout << "Perdas = " << fo << endl;
			if (m_verbose) {
				cout << "Perdas = " << fo << endl;
			}
			m_foMod(k) = m_calc->Fitness(graph, m_w);
			m_perdas(k) = fo;
			//break;
			// Passo 5 - Teste de convergência:
			int nElem = (int) dLdu.n_rows;
			vec gradRed = zeros<vec>(nElem);
			gradRed = dLdu;

			int size = (int) graph->GetNumBus();
			for (int i = 0; i < size; i++) {
				Bus* bus = graph->GetBus(i + 1);
				if (bus->GetType() != Bus::LOAD) {
					if (bus->GetVCalc() - m_erro <= bus->GetBus().m_vmin) {
						if (gradRed(bus->GetOrdG()) >= 0) {
							gradRed(bus->GetOrdG()) = 0;
						}
					} else if (bus->GetVCalc() + m_erro >= bus->GetBus().m_vmax) {
						if (gradRed(bus->GetOrdG()) <= 0) {
							gradRed(bus->GetOrdG()) = 0;
						}
					}
				}
			}

			if (m_verbose) {
				cout << "Teste de convergência: " << endl;
			}
			if (norm(gradRed) <= m_erro) {
				conv = 1;
			} else {
				conv = 0;
			}

			if (conv == 1) {
				flag = 1;
			} else if (k != m_maxIter) {
				// Passo 6 - Ajustar as variáveis de controle u, fazer
				// k = k+1 e voltar para o passo 2:
				k = k + 1;
				m_w = m_w * m_ro;
				cout << m_w << endl;

				double c = m_co;
				Data_t bu = m_calc->Busca(lf, graph, m_w, c, m_erro, dLdu);
				break;
				if (bu.m_conv == 1) {
					// Ajuste das magnitudes de tensão das barras de geração:
					// V = ajuste_controles(npv, V, V_min, V_max, busG, c, dLdu, 1);
					for (int i = 0; i < size; i++) {
						Bus* bus = graph->GetBus(i + 1);
						if (bus->GetType() != Bus::GENERATION) {
							bus->SetVCalc(bus->GetVCalc() + c * dLdu(i));
						}
					}
					// Correção dos limites violados:
					for (int i = 0; i < size; i++) {
						Bus* bus = graph->GetBus(i + 1);
						if (bus->GetType() != Bus::GENERATION
								&& bus->GetVCalc() < bus->GetBus().m_vmin) {
							bus->SetVCalc(bus->GetBus().m_vmin);
							if (m_verbose) {
								printf("\nV%d < V_min! -> V%d = V_min",
										bus->GetVCalc(), bus->GetVCalc());
							}
						} else if (bus->GetType() != Bus::GENERATION
								&& bus->GetVCalc() > bus->GetBus().m_vmax) {
							bus->SetVCalc(bus->GetBus().m_vmax);
							if (m_verbose) {
								printf("\nV%d > V_max! -> V%d = V_max",
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

					//  Armazenamento das variáveis de controle a cada iteração:
					for (int i = 0; size; i++) {
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
}

// Relatório de saída:

if (m_verbose) {
	if (flag == 1) {
		printf(
				"\n\nO método do Gradiente Reduzido convergiu em %d iterações com E = %.0g!",
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
		printf("\n\nResumo do Relatório de Convergência:");
		/*printf("\nIteração   V1        V2        FO Mod    Perdas");
		 printf("[it]       [p.u.]    [p.u.]    [p.u.]    [p.u.]");
		 for i=1:k+1,
		 disp(sprintf('%3d %13.4f %9.4f %9.4f %9.4f', i-1, ctrls(i, 1), ctrls(i, 2), fo_mod(i), perdas(i));
		 end*/

		printf("\n\nRelatório de Saída:");
		Report* rp = new Report();
		rp->Output(graph);
	}
}
}

}
