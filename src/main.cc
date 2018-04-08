#include "iostream"
#include "armadillo"
#include "fpo/fpo.h"
#include "pso/pso.h"
#include "pso/sphere.h"

using namespace arma;
using namespace load;
using namespace fpo;
using namespace pso;

int main(int argc, char **argv) {
	/*std::cout << "Rodando o Fluxo" << std::endl;
	Fpo* f = new Fpo;
	f->Execute("/home/pc-thiago/eclipse-workspace/LoadFlow/src/data/3bus.txt");
	delete f;*/
	Problem* p = new Sphere();
	p->SetNVar(5);
	p->SetVarMin(-10 * ones<vec>(5));
	p->SetVarMax(10 * ones<vec>(5));
	p->SetVarSize(zeros<vec>(5));
	p->SetType(MIN);
	p->SetVarSize(zeros<vec>(5));

	Pso* pso = new Pso();
	pso->SetProblem(p);
	pso->Execute();
	delete pso;
	/*cout << "Armadillo version: " << arma_version::as_string() << endl;

	mat A;

	A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << endr
			<< 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << endr
			<< 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << endr
			<< 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << endr
			<< 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << endr;

	A.print("A =");

	// determinant
	cout << "det(A) = " << det(A) << endl;

	// inverse
	cout << "inv(A) = " << endl << inv(A) << endl;*/
}
