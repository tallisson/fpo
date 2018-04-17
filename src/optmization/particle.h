#ifndef OPT_PARTICLE_H_
#define OPT_PARTICLE_H_

#include <armadillo>
#include <string>

using namespace arma;
using namespace std;

namespace opt {

class Particle {
public:
	Particle();
	virtual ~Particle();

	void SetPosition(vec position);
	vec GetPosition(void) const;

	void SetVelocity(vec velocity);
	vec GetVelocity(void) const;

	void SetCost(double cost);
	double GetCost(void) const;

	void SetBestPos(vec bestPos);
	vec GetBestPos(void) const;

	void SetBestCost(double bestCost);
	double GetBestCost(void) const;

	vec GetMaxVel(void) const;
	void SetMaxVel(vec maxVel);

	vec GetMinVel(void) const;
	void SetMinVel(vec minVel);

	void Print(void);
private:
	static int m_qtdParticles;
	int m_id;
	vec m_position;
	vec m_velocity;
	double m_cost;
	vec m_bestPos;
	double m_bestCost;
	vec m_maxVel;
	vec m_minVel;

};

} /* namespace opt */

#endif /* OPT_PARTICLE_H_ */
