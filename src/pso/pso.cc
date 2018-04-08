#include "pso.h"

#include <math.h>
#include <armadillo>
#include <iostream>
#include <iostream>
#include <fstream>

using namespace arma;
using namespace std;

namespace pso {

double Pso::INF = 10000000;
double Pso::INF_NEG = -10000000;

Pso::Pso() {
	// Constraints Coefficients
	m_kappa = 1;
	m_phi1 = 2.05;
	m_phi2 = 2.05;
	m_phi = m_phi1+m_phi2;
	m_chi = 2 * m_kappa/ abs(2-m_phi-sqrt(pow(m_phi, 2)-4*m_phi));

	m_maxIt = 1000;
	m_nPop = 50;
	m_w = 1;
	m_wdamp = 0.99;
	m_c1 = 2;
	//m_c1 = chi*phi1;
	m_c2 = 2;
	//c2 = chi*phi2;
	m_globalBest = new Particle();
	m_globalBest->SetCost(INF);
	m_problem = 0;
}

Pso::~Pso() {
	m_particles.clear();
	delete m_globalBest;
	m_globalBest = 0;
	delete m_problem;
	m_problem = 0;
}

double Pso::GetKappa(void) const {
	return m_kappa;
}

void Pso::SetKappa(double kappa) {
	m_kappa = kappa;
}

double Pso::GetPhi1(void) const {
	return m_phi;
}

void Pso::SetPhi1(double phi1) {
	m_phi1 = phi1;
}

double Pso::GetPhi2(void) const {
	return m_phi2;
}

void Pso::SetPhi2(double phi2) {
	m_phi2 = phi2;
}

double Pso::GetPhi(void) const {
	return m_phi;
}

void Pso::SetPhi(double phi) {
	m_phi = phi;
}

int Pso::GetMaxIt(void) const {
	return m_maxIt;
}

void Pso::SetMaxIt(int maxIt) {
	m_maxIt = maxIt;
}

double Pso::GetW(void) const {
	return m_w;
}

void Pso::SetW(double w) {
	m_w = w;
}

double Pso::GetWDamp(void) const {
	return m_wdamp;
}

void Pso::SetWDamp(double wdamp) {
	m_wdamp = wdamp;
}

double Pso::GetC1(void) const {
	return m_c1;
}

void Pso::SetC1(double c1) {
	m_c1 = c1;
}

double Pso::GetC2(void) const {
	return m_c2;
}

void Pso::SetC2(double c2) {
	m_c2 = c2;
}

Problem* Pso::GetProblem(void) const {
	return m_problem;
}

void Pso::SetProblem(Problem* problem) {
	if(problem->GetType() == MAX) {
		m_globalBest->SetCost(INF_NEG);
	}
	m_problem = problem;
}

void Pso::CreatePopulation(void) {
	for(int i = 0; i < m_nPop; i++) {
		Particle* p = new Particle();
		m_particles.push_back(p);
	}
	InitPop();
	m_bestCosts = zeros<vec>(m_maxIt);
}

// Initialize Population Members
void Pso::InitPop(void) {
	if(m_problem == 0) {
		cout << "Especificar problema para o pso (132->Pso.cc)" << endl;
		return;
	}
	for(int n = 0; n < m_nPop; n++) {
		Particle* particle = m_particles.at(n);
		// Generation Random Solution
		vec position = zeros<vec>(m_problem->GetNVar());
		vec max = m_problem->GetVarMax();
		vec min = m_problem->GetVarMin();
		arma_rng::set_seed_random();
		for(int i = 0; i < m_problem->GetNVar(); i++) {
			double x = (max.at(i) - min.at(i)) * randu<double>() + min.at(i);
			position(i) = x;
		}
		particle->SetPosition(position);

		// Initialize Velocity
		particle->SetVelocity(zeros<vec>(m_problem->GetNVar()));

		// Evaluation
		double cost = m_problem->CostFunction(position);
		particle->SetCost(cost);

		// Update the Personal Best
		particle->SetBestPos(position);
		particle->SetBestCost(cost);

		// Update Global Best
		if (m_problem->Compare(particle->GetBestCost(), m_globalBest->GetCost())) {
			m_globalBest->SetPosition(particle->GetBestPos());
			m_globalBest->SetCost(particle->GetBestCost());
		}
		particle->Print();
	}
	m_globalBest->Print();
}

// Main Loop of PSO
void Pso::MainLoop(void) {
	arma_rng::set_seed_random();
	for(int it = 0; it < m_maxIt; it++) {

	  for(int i = 0; i < m_nPop; i++) {
	    // Update Velocity
		Particle* particle = m_particles.at(i);

		// Inertial Component
		vec inertialComponent = zeros<vec>(m_problem->GetNVar());
		inertialComponent = m_w * particle->GetVelocity();

		// Calculate Cognitive Component
		vec c1AndVarSize = m_c1 * randu<vec>(m_problem->GetNVar());
		vec diffPosAndBest = particle->GetBestPos() - particle->GetPosition();
		vec cognitiveComponent = zeros<vec>(m_problem->GetNVar());
		cognitiveComponent = c1AndVarSize % diffPosAndBest;

		// Calculate Social Component
		vec c2AndVarSize = m_c2 * randu<vec>(m_problem->GetNVar());
		vec diffPosAndGlobal = m_globalBest->GetPosition() - particle->GetPosition();
		vec socialComponent = zeros<vec>(m_problem->GetNVar());
		socialComponent = c2AndVarSize % diffPosAndGlobal;

		vec newVel = zeros<vec>(m_problem->GetNVar());
		newVel = inertialComponent + cognitiveComponent + socialComponent;

		particle->SetVelocity(newVel);
	    /* Apply velocity Limits
	       particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
	       particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
	     */

	    // Update Position
	    particle->SetPosition(particle->GetPosition() + newVel);

	    /* Apply Lower and Upper Bound Limits
	       particle(i).Position = max(particle(i).Position, VarMin);
	       particle(i).Position = min(particle(i).Position, VarMax);
	     */

	    // Update Cost
	    double cost = m_problem->CostFunction(particle->GetPosition());
	    particle->SetCost(cost);

	    if(m_problem->Compare(particle->GetCost(), particle->GetBestCost())) {
	    	particle->SetBestPos(particle->GetPosition());
            particle->SetBestCost(particle->GetCost());

			// Update Global Best
			if(particle->GetBestCost() < m_globalBest->GetCost()) {
			  m_globalBest->SetPosition(particle->GetBestPos());
			  m_globalBest->SetCost(particle->GetBestCost());
			}
	    }

	    /*if(particle->GetCost() < particle->GetBestCost()) {

	      particle->SetBestPos(particle->GetPosition());
	      particle->SetBestCost(particle->GetCost());

	      // Update Global Best
	      if(particle->GetBestCost() < m_globalBest->GetCost()) {
	        m_globalBest->SetPosition(particle->GetBestPos());
	        m_globalBest->SetCost(particle->GetBestCost());
	      }

	    }*/
	  }
	  // Store the Best Cost Value
	  double gCost =
	  m_bestCosts.at(it) = m_globalBest->GetCost();

	  // Display Iteration Information
	  cout << "[ Iteration " << (it+1) << " : Best Cost = " << gCost << " ]" << endl;

	  // Damping Inertial Coefficient
	  m_w *= m_wdamp;
	}
	SaveResults();
}

void Pso::Execute(void) {
	CreatePopulation();
	MainLoop();
}

void Pso::SaveResults(void) {
  ofstream file;
  file.open ("data.txt");
  int size = (int)m_bestCosts.size();
  for(int i = 0; i < size; i++) {
	  file << (i+1) << "\t" << m_bestCosts.at(i) << endl;
  }
  file << "EOF";
  file.close();
}

} /* namespace pso */
