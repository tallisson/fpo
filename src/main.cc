#include "iostream"
#include <armadillo>
#include "fpo/fpo.h"
#include "pso/pso.h"
#include "pso/sphere.h"
#include "optmization/fpo-pso.h"
#include "fpo/fpo-reactive.h"

#include <set>

using namespace arma;
using namespace load;
using namespace fpo;
using namespace pso;
using namespace std;
using namespace opt;

int main(int argc, char **argv) {
	std::cout << "Rodando o Fluxo" << std::endl;

	/*FpoPso* f = new FpoPso();
	f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/data/6bus.txt");
	delete f;
	f = 0;*/

	FpoReactive* f = new FpoReactive();
	f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/data/6bus.txt");
	delete f;
	f = 0;
	return 1;
}
