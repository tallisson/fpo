#include "fpo-pso.h"
#include "fpo-problem.h"
#include "../load/calc.h"

using namespace load;
using namespace opt;
using namespace model;

namespace opt {

FpoPso::FpoPso() {
	// TODO Auto-generated constructor stub
	m_lf = new LoadFlow();
	m_pso = 0;
	m_maxIter = 100;
	m_w = 85;
	m_damp = 0.99;
	m_ro = 1;
	m_co = 0.00124;
	m_verbose = false;
	m_pso = new Pso();
	m_foMod = zeros(m_maxIter + 1);
	m_perdas = zeros(m_maxIter + 1);
}

void FpoPso::InitPso(string cdf, double cost) {
	//FpoProblem* p = new FpoProblem(cdf, m_w);
	FpoProblem* p = new FpoProblem(m_lf, m_w);
	// Definindo o número de barras de Geração
	int numG = m_lf->GetGraph()->GetNumPV() + 1;
	p->SetNVar(numG);

	// limites máximo
	vec varMax = zeros<vec>(numG);
	vec varMin = zeros<vec>(numG);
	int numB = m_lf->GetGraph()->GetNumBus();

	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);
		if (bus->GetType() == Bus::GENERATION || bus->GetType() == Bus::SLACK) {
			varMax.at(bus->GetOrdG()) = 1.05;/*bus->GetBus().m_vmax;*/
			varMin.at(bus->GetOrdG()) = 0.95;/*bus->GetBus().m_vmin;*/
		}
	}

	p->SetVarMin(varMin);
	p->SetVarMax(varMax);
	p->SetVarSize(zeros<vec>(numG));
	p->SetType(MIN);

	m_pso->SetProblem(p);
	vec pos = zeros<vec>(numG);
	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);
		if (bus->GetType() == Bus::GENERATION || bus->GetType() == Bus::SLACK) {
			pos.at(bus->GetOrdG()) = bus->GetVCalc();
		}
	}
	m_pso->ParticleBase(pos, cost);
}

FpoPso::~FpoPso() {
	// TODO Auto-generated destructor stub
	m_x.clear();
	m_foMod.clear();
	m_perdas.clear();
	delete m_lf;
	m_lf = 0;
	delete m_pso;
	m_pso = 0;
}

void FpoPso::Execute(std::string cdf) {
	m_lf->Prepare(cdf);
	Calc* calc = new Calc();

	Graph* graph = m_lf->GetGraph();
	for (int i = 0; i < graph->GetNumBus(); i++) {
		Bus* bus = graph->GetBus(i + 1);
		if (bus->GetType() != Bus::LOAD
				&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			bus->SetVCalc(1.0);
		}
	}

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

		int conv = m_lf->Execute();
		if (m_verbose) {
			for (int i = 1; i <= size; i++) {
				Bus* bus = graph->GetBus(i);
				printf("%.12f\n", bus->GetVCalc());
			}
			printf("\n");
		}

		if (conv == 1) {
			//problem->SetW(0);
			double fo = calc->Fitness(m_lf->GetGraph(), 0);

			if (m_verbose) {
				cout << "Perdas = " << fo << endl;
			}
			//problem->SetW(m_w);
			m_foMod(k) = calc->Fitness(m_lf->GetGraph(), m_w);
			cout << m_foMod(k) << endl;
			InitPso(cdf, m_foMod(k));

			m_perdas(k) = fo;
			//cout << "Perdas = " << fo;
			if (k++ != m_maxIter) {
				//m_w *= m_damp;
				int f = m_pso->Execute();
				cout << "F = " << f << endl;
				Particle* p = m_pso->GetBest();
				cout << "Best Particle = " << endl;
				p->Print();

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

				for (int i = 0; i < m_lf->GetGraph()->GetNumBus(); i++) {
					Bus* bus = graph->GetBus(i + 1);
					if (bus->GetType() != Bus::GENERATION) {
						bus->AddControl(bus->GetVCalc());
					}
				}

				if (m_verbose) {
					printf("\n\n-----------------------------------------");
					printf("\nIteração %d:", k);
				}

				if (f != 1) {
					flag = 1;
				}
			}
		} else {
			flag = 0;
		}
	}

	// Relatório de saída:
	m_verbose = true;
	if (m_verbose) {
		if (flag == 1) {
			printf(
					"\n\nO método do Gradiente Reduzido convergiu em %d iterações",
					k/*, m_erro*/);
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
			Report(k);
		}
	}

}

void FpoPso::UpdateV(vec v) {
	int numB = m_lf->GetGraph()->GetNumBus();
	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);
		if (bus->GetType() != Bus::LOAD
				&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			bus->SetVCalc(v.at(bus->GetOrdG()));
		}
	}
}

void FpoPso::Report(int numIt) {
	printf("\n\nResumo do Relatório de Convergência:");
	printf("\n\nRelatório de Saída:");
	printf("\nIteração   V1       V2        FO Mod    Perdas\n");
	printf("[it]       [p.u.]    [p.u.]    [p.u.]    [p.u.]\n");
	for (int i = 0; i <= numIt; i++) {
	  printf("%3d %13.4f %9.4f %9.4f %9.4f\n", i, 0.0, 0.0, m_foMod.at(i),m_perdas.at(i));
	}
	printf("\n\nRelatório de Saída:\n");
	int numB = (int) m_lf->GetGraph()->GetNumBus();
	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);
		bus->Print();
	}

	// Geração de potência ativa:
	Bus* slack = m_lf->GetGraph()->GetSlack();
	slack->CalcPG();

	// Geração de potência reativa:
	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);

		if (bus->GetType() != Bus::LOAD
				&& bus->GetType() != Bus::LOSS_CONTROL_REACT) {
			bus->CalcQg();
		}
	}

	// Impressão na tela do relatório de saída:

	printf(
			"\nBarra   V         ang       ang       Pg        Qg        Pc        Qc\n");
	printf(
			"[k]     [p.u.]    [graus]   [rad]     [p.u.]    [p.u.]    [p.u.]    [p.u.]\n");
	for (int i = 0; i < numB; i++) {
		Bus* bus = m_lf->GetGraph()->GetBus(i + 1);
		printf("%2d %11.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
				bus->GetBus().m_nin, bus->GetVCalc(),
				bus->GetACalc() * (180 / M_PI), bus->GetACalc(),
				bus->GetPCalc(), bus->GetQCalc(), bus->GetBus().m_pc,
				bus->GetBus().m_qc);
	}

	// Cálculo dos fluxos nas linhas de transmissão e perdas ativas:
	DataLoss_t data = m_lf->CalcLosses();
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


void FpoPso::SetVerbose(bool verbose) {
	m_verbose = verbose;
}

void FpoPso::SetW(double w) {
	m_w = w;
}

} /* namespace load */
