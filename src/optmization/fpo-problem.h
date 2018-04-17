#ifndef OPT_FPO_PROBLEM_H_
#define OPT_FPO_PROBLEM_H_

#include "fpo-problem.h"
#include "../load/load-flow.h"
#include "particle.h"

#include <armadillo>
#include <string>

using namespace arma;
using namespace load;
using namespace model;
using namespace std;

namespace opt {

struct vLim {
	vec m_max;
	vec m_min;

	vLim(vec min, vec max) {
		m_max = max;
		m_min = min;
	}
};
typedef vLim VLim;

enum topt {
	MIN = 0, MAX = 1
};
typedef topt TypeOpt;

struct aresta {
	aresta(int k, int m) :
			m_k(k), m_m(m) {
	}
	int m_k;
	int m_m;

	bool equals(aresta rhs) {
		if (m_k == rhs.m_k && m_m == rhs.m_m) {
			return true;
		}
		if (m_k == rhs.m_m && m_m == rhs.m_k) {
			return true;
		}
		return false;
	}
};
typedef aresta Aresta_t;

class FpoProblem {
public:
	FpoProblem();
	FpoProblem(string cdf, double w);
	FpoProblem(LoadFlow* lf, double w);

	virtual ~FpoProblem();
	double CostFunction(Particle* p);

	double GetNVar() const;
	void SetNVar(double nVar);

	vec GetVarMax() const;
	void SetVarMax(vec varMax);

	vec GetVarMin() const;
	void SetVarMin(vec varMin);

	vec GetVarSize() const;
	void SetVarSize(vec varSize);

	vec GetMaxVelocity(void) const;
	vec GetMinVelocity(void) const;

	VLim CalcVelLimits(void);

	TypeOpt GetType(void) const;
	void SetType(TypeOpt type);

	bool Compare(double value1, double value);

	void SetW(double w);

	void SetCdf(string cdf);
	void Print(void);
private:
	double m_w;
	LoadFlow* m_lf;
	string m_cdf;

	double m_nVar; // Number of Unknown (Decision) Variables
	vec m_varSize; // Matrix Size of Decision Variables [1 nVar]
	vec m_varMin; // Lower Bound of Decision Variables
	vec m_varMax; // Upper Bound of Decision Variables

	TypeOpt m_type;

	vec m_maxVelocity;
	vec m_minVelocity;
};

} /* namespace opt */

#endif /* OPT_FPO_PROBLEM_H_ */
