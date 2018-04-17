#include "iostream"
#include "armadillo"
#include "fpo/fpo.h"
#include "pso/pso.h"
#include "pso/sphere.h"
#include "optmization/fpo-pso.h"
#include <set>

using namespace arma;
using namespace load;
using namespace fpo;
using namespace pso;
using namespace std;
using namespace opt;

int main(int argc, char **argv) {
	std::cout << "Rodando o Fluxo" << std::endl;
	//FpoPso* f = new FpoPso;
	//f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/data/6bus.txt");
	//Fpo* f = new Fpo;
	//f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/data/6bus.txt");
	//delete f;
	/*Problem* p = new Sphere();
	p->SetNVar(5);
	p->SetVarMin(-10 * ones<vec>(5));
	p->SetVarMax(10 * ones<vec>(5));
	p->SetVarSize(zeros<vec>(5));
	p->SetType(MIN);
	p->SetVarSize(zeros<vec>(5));

	Pso* pso = new Pso();
	pso->SetProblem(p);
	pso->Execute();
	delete pso;*/
	FpoPso* f = new FpoPso();
	f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/data/6bus.txt");
	return 1;
}
