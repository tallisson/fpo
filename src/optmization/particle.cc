#include "particle.h"
#include <iostream>

namespace opt {

int Particle::m_qtdParticles = 0;

Particle::Particle() {
	m_cost = 0;
	m_bestCost = 0;
	m_id = 0;
	m_id = m_qtdParticles++;
}

Particle::~Particle() {
}

void Particle::SetPosition(vec position) {
	m_position = position;
}

vec Particle::GetPosition(void) const {
	return m_position;
}

void Particle::SetVelocity(vec velocity) {
	m_velocity = velocity;
}

vec Particle::GetVelocity(void) const {
	return m_velocity;
}

void Particle::SetCost(double cost) {
	m_cost = cost;
}

double Particle::GetCost(void) const {
	return m_cost;
}

void Particle::SetBestPos(vec bestPos) {
	m_bestPos = bestPos;
}

vec Particle::GetBestPos(void) const {
	return m_bestPos;
}

double Particle::GetBestCost(void) const {
	return m_bestCost;
}

void Particle::SetBestCost(double bestCost) {
	m_bestCost = bestCost;
}

vec Particle::GetMaxVel(void) const {
	return m_maxVel;
}

void Particle::SetMaxVel(vec maxVel) {
	m_maxVel = maxVel;
}

vec Particle::GetMinVel(void) const {
	return m_minVel;
}
void Particle::SetMinVel(vec minVel) {
	m_minVel = minVel;
}

void Particle::Print(void) {
	cout << "Particle " << m_id << ":" << endl
		 << "Position = " << endl << m_position
		 << "Cost = " << m_cost << endl
		 << "+++++++++++++++++++++++++++++++++++++" << endl;
}

} /* namespace pso */
