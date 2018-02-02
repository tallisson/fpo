/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 */

#include "branch.h"

#include <iostream>
#include <math.h>

namespace model {

Branch::Branch() {
}

Branch::~Branch() {
}

void Branch::SetBranch(DBranch_t branch) {
	m_branch = branch;
}

DBranch_t Branch::GetBranch() const {
	return m_branch;
}

double Branch::CalcPkmL(double vK, double vM, double aK, double aM) {
	// s.branch.g(i)*s.bus.v(k)^2 - s.bus.v(k)*s.bus.v(m)*(s.branch.g(i)*cos(akm)+s.branch.b(i)*sin(akm));
	double theta_km = aK - aM;
	m_p_km_L = GetBranch().m_g * pow(vK, 2)
			- vK * vM
					* (GetBranch().m_g * cos(theta_km)
							+ GetBranch().m_b * sin(theta_km));

	return m_p_km_L;
}

double Branch::CalcPmkL(double vK, double vM, double aK, double aM) {
	// s.branch.g(i)*s.bus.v(m)^2 - s.bus.v(k)*s.bus.v(m)*(s.branch.g(i)*cos(akm)-s.branch.b(i)*sin(akm));
	double theta_km = aK - aM;
	m_p_mk_L = GetBranch().m_g * pow(vM, 2)
			- vK * vM
					* (GetBranch().m_g * cos(theta_km)
							- GetBranch().m_b * sin(theta_km));

	return m_p_mk_L;
}

double Branch::CalcQkmL(double vK, double vM, double aK, double aM) {
	// -(s.branch.b(i)+s.branch.bsh(i))*s.bus.v(k)^2 + s.bus.v(k)*s.bus.v(m)*(s.branch.b(i)*cos(akm)-s.branch.g(i)*sin(akm))
	double theta_km = aK - aM;
	m_q_km_L = -(GetBranch().m_b + GetBranch().m_bsh) * pow(vK, 2)
			+ vK * vM
					* (GetBranch().m_b * cos(theta_km)
							- GetBranch().m_g * sin(theta_km));

	return m_q_km_L;
}

double Branch::CalcQmkL(double vK, double vM, double aK, double aM) {
	// -(s.branch.b(i)+s.branch.bsh(i))*s.bus.v(m)^2 + s.bus.v(k)*s.bus.v(m)*(s.branch.b(i)*cos(akm)+s.branch.g(i)*sin(akm));
	double theta_km = aK - aM;
	m_q_mk_L = -(GetBranch().m_b + GetBranch().m_bsh) * pow(vM, 2)
			+ vK * vM
					* (GetBranch().m_b * cos(theta_km)
							+ GetBranch().m_g * sin(theta_km));

	return m_q_mk_L;
}

double Branch::CalcL(double vK, double vM, double aK, double aM) {
	// s.branch.g(i)*(s.bus.v(k)^2 + s.bus.v(m)^2 - 2*s.bus.v(k)*s.bus.v(m)*cos(akm));
	double theta_km = aK - aM;
	m_l = GetBranch().m_g
			* (pow(vK, 2) + pow(vM, 2)
					- 2 * vK * vM * cos(theta_km));

	return m_l;
}

void Branch::Print(void) {
	std::cout << "Branch ( " << m_branch.m_ni << ", " << m_branch.m_nf
			<< ") => " << m_branch.m_b << "\t" << m_branch.m_bsh << "\t"
			<< m_branch.m_g << "\t" << m_branch.m_tap << "\t" << m_branch.m_def
			<< std::endl;
}

}
