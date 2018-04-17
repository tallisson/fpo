#ifndef OPT_PSO_H_
#define OPT_PSO_H_

#include "particle.h"
#include "fpo-problem.h"
#include <vector>

using namespace std;

namespace opt {

class Pso {
public:
	Pso();
	virtual ~Pso();

	double GetKappa(void) const;
	void SetKappa(double kappa);

	double GetPhi1(void) const;
	void SetPhi1(double phi1);

	double GetPhi2(void) const;
	void SetPhi2(double phi2);

	double GetPhi(void) const;
	void SetPhi(double phi);

	int GetMaxIt(void) const;
	void SetMaxIt(int maxIt);

	double GetW(void) const;
	void SetW(double w);

	double GetWDamp(void) const;
	void SetWDamp(double wdamp);

	double GetC1(void) const;
	void SetC1(double c1);

	double GetC2(void) const;
	void SetC2(double c2);

	Particle* GetBest(void) const;

	FpoProblem* GetProblem(void) const;
	void SetProblem(FpoProblem* problem);

	// Create Population
	void CreatePopulation(void);

	void MainLoop(void);

	int Execute(void);

	void SaveResults(void);

	void ParticleBase(vec position, double cost);

	static double INF;
	static double INF_NEG;
private:
	double m_kappa;
	double m_phi1;
	double m_phi2;
	double m_phi;
	double m_chi;

	int m_maxIt; // Maximium Number of Iterations
	double m_nPop; // Population Size (Swarm Size)
	double m_w; // Inertial Coefficient
	double m_wdamp; // Damping Ratio of Inertia Coefficient
	double m_c1; // Personal Acceleration Coefficient
	double m_c2; // Social Acceleration Coefficient
	int m_execNum;

	// Array to Hold Best Cost Value on Each Iteration
	vec m_bestCosts;
	vec m_allCosts;

	Particle* m_globalBest;
	vector<Particle* > m_particles;
	FpoProblem* m_problem;

	void InitPop(void);
};

} /* namespace opt */

#endif /* OPT_PSO_H_ */
