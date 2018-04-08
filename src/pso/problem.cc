#include "problem.h"

namespace pso {

Problem::Problem() {
	m_nVar = 0;
	m_type = MIN;
}

Problem::~Problem() {
	// TODO Auto-generated destructor stub
}

double Problem::GetNVar() const {
	return m_nVar;
}

void Problem::SetNVar(double nVar) {
	m_nVar = nVar;
}

vec Problem::GetVarMax() const {
	return m_varMax;
}

void Problem::SetVarMax(vec varMax) {
	m_varMax = varMax;
}

vec Problem::GetVarMin() const {
	return m_varMin;
}

void Problem::SetVarMin(vec varMin) {
	m_varMin = varMin;
}

vec Problem::GetVarSize() const {
	return m_varSize;
}

void Problem::SetVarSize(vec varSize) {
	m_varSize = varSize;
}

vec Problem::GetMaxVelocity(void) const {
	return m_maxVelocity;
}

vec Problem::GetMinVelocity(void) const {
	return m_minVelocity;
}

TypeOpt Problem::GetType(void) const {
	return m_type;
}

void Problem::SetType(TypeOpt type) {
	m_type = type;
}

VLim Problem::CalcVelLimits(void) {
	m_maxVelocity = 0.2*(m_varMax - m_varMin);
	m_minVelocity = -m_maxVelocity;
	VLim lim(m_minVelocity, m_maxVelocity);
	return lim;
}

bool Problem::Compare(double value1, double value2) {
	switch (m_type) {
		case MIN:
			return value1 < value2;
			break;
		case MAX:
			return value1 > value2;
			break;
		default:
			return value1 > value2;
			break;
	}
}

} /* namespace pso */
