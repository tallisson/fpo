#ifndef OPT_FPO_PSO_H_
#define OPT_FPO_PSO_H_

#include "../load/load-flow.h"
#include "fpo-pso.h"
#include "pso.h"
#include <armadillo>
#include <string>

using namespace arma;
using namespace load;
using namespace std;

namespace opt {

class FpoPso{
public:
	FpoPso();
	virtual ~FpoPso();

	void SetVerbose(bool verbose);
	void SetW(double w);

	void Execute(std::string cdf);

	void Report(int numIt);
private:
	bool m_verbose;
	double m_w;
	double m_damp;
	int m_co;
	int m_ro;
	int m_maxIter;
	vec m_foMod;
	vec m_perdas;
	vec m_x;
	LoadFlow* m_lf;
	Pso* m_pso;

	void InitPso(string cdf, double cost);
	void UpdateV(vec v);
};

} /* namespace fpo */

#endif /* OPT_FPO_PSO_H_ */
