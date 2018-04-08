#ifndef PROBLEM_H_
#define PROBLEM_H_

#include <armadillo>

using namespace arma;

namespace pso {

struct vLim {
	vec m_max;
	vec m_min;

	vLim(vec min, vec max) {
		m_max = max;
		m_min = min;
	}
};
typedef vLim VLim;

enum topt{
	MIN = 0, MAX = 1
};
typedef topt TypeOpt;

class Problem {
public:
	Problem();
	virtual ~Problem();
	virtual double CostFunction(vec x) = 0;

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
private:
	double m_nVar; // Number of Unknown (Decision) Variables
	vec m_varSize; // Matrix Size of Decision Variables [1 nVar]
	vec m_varMin; // Lower Bound of Decision Variables
	vec m_varMax; // Upper Bound of Decision Variables

	TypeOpt m_type;

	vec m_maxVelocity;
	vec m_minVelocity;
};

} /* namespace pso */

#endif /* PROBLEM_H_ */
