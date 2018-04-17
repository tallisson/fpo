/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "bus.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <sstream>

namespace model {

Bus::Bus() {
	m_tap = IMP;
	m_crt = 0.0;
	m_dsv = 0;
	m_vCalc = 0.0;
	m_aCalc = 0.0;
	m_isControlled = false;
	m_status = NORMAL;
	m_type = NONE;
	m_pgCalc = 0;
	m_qgCalc = 0;
	m_vCalc = 0;
	m_aCalc = 0;
	m_ordG = -1;
	m_numControl = 0;
	m_vPreBusca = 0;
	m_aPreBusca = 0;
}

Bus::~Bus() {
}

void Bus::SetBus(DBus_t bus) {
	m_bus = bus;
}

DBus_t Bus::GetBus() const {
	return m_bus;
}

void Bus::SetType(Bus::Type type) {
	m_type = type;
}

Bus::Type Bus::GetType(void) {
	return m_type;
}

void Bus::AddBranch(Branch* branch, Bus* neighbor) {
	m_branches.push_back(branch);
	m_neighbors.push_back(neighbor);
}

std::vector<Branch*> Bus::GetBranches() const {
	return m_branches;
}

double Bus::CalcPg(void) {
	double pgK = 0;
	int size = m_branches.size();
	for (int i = 0; i < size; i++) {
		Branch* branch = m_branches.at(i);
		Bus* busM = m_neighbors.at(i);
		DBranch_t dataBranch = branch->GetBranch();

		double aM = busM->GetACalc();
		double vM = busM->GetVCalc();
		double theta_km = m_aCalc - aM;
		double vK = m_vCalc;

		if (m_bus.m_nin == dataBranch.m_nf) {
			pgK += (dataBranch.m_g * (1 / pow(dataBranch.m_tap, 2)) * pow(vK, 2)
					- (1 / dataBranch.m_tap) * vK * vM * (dataBranch.m_g * cos(theta_km - dataBranch.m_def)
					+ dataBranch.m_b * sin(theta_km - dataBranch.m_def)));
		} else {
			pgK += (dataBranch.m_g * pow(vK, 2) - (1 / dataBranch.m_tap) * vK * vM
				   * (dataBranch.m_g * cos(theta_km + dataBranch.m_def) + dataBranch.m_b
				   * sin(theta_km + dataBranch.m_def)));
		}
	}
	pgK += -GetBus().m_gsh * pow(m_vCalc, 2) + GetBus().m_pc;
	m_pgCalc = pgK;

	return m_pgCalc;
}

double Bus::CalcQg(void) {
	double qgK = 0;
	int size = m_branches.size();
	for (int i = 0; i < size; i++) {
		Branch* branch = m_branches.at(i);
		Bus* busM = m_neighbors.at(i);
		DBranch_t dataBranch = branch->GetBranch();

		double aM = busM->GetACalc();
		double vM = busM->GetVCalc();
		/*busM->SetACalc(aM);
		busM->SetVCalc(vM);*/
		double theta_km = m_aCalc - aM;
		double vK = m_vCalc;

		if (m_bus.m_nin == dataBranch.m_ni) {
			qgK += (-(dataBranch.m_b * (1 / pow(dataBranch.m_tap, 2)) + dataBranch.m_bsh) * pow(vK, 2)
				   + (1 / dataBranch.m_tap) * vK * vM * (dataBranch.m_b * cos(theta_km - dataBranch.m_def)
				   - dataBranch.m_g * sin(theta_km - dataBranch.m_def)));
		} else {
			qgK += (-dataBranch.m_b * pow(vK, 2) + (1 / dataBranch.m_tap) * vK * vM
				   * (dataBranch.m_b * cos(theta_km + dataBranch.m_def) - dataBranch.m_g
				   * sin(theta_km + dataBranch.m_def)));
		}
	}
	qgK += -GetBus().m_bsh * pow(m_vCalc, 2) + GetBus().m_qc;
	m_qgCalc = qgK;

	return m_qgCalc;
}

std::vector<Bus*> Bus::GetNeighbors(void) const {
	return m_neighbors;
}

void Bus::SetTap(Bus::TapType tap) {
	m_tap = tap;
}

Bus::TapType Bus::GetTap(void) const {
	return m_tap;
}

void Bus::Print(void) {
	std::cout.precision(11);

	std::cout << "Bus = " << m_bus.m_nin << "\t" << "V = " << m_vCalc << "\t"
			<< "ang = " << m_aCalc << std::endl;
}

void Bus::SetCrt(double crt) {
	m_crt = crt;
}

double Bus::GetCrt(void) const {
	return m_crt;
}

double Bus::CalcDsv(void) {
	m_dsv = 0;
	if (m_vCalc < MIN_VOLTAGE_ONS) {
		m_dsv = Bus::MIN_VOLTAGE_ONS - m_vCalc;
		m_status = MIN_VOLTAGE_VIOLATION;
		m_isControlled = true;
	}

	if (m_vCalc > Bus::MAX_VOLTAGE_ONS) {
		m_dsv = m_vCalc - Bus::MAX_VOLTAGE_ONS;
		m_status = MAX_VOLTAGE_VIOLATION;
		m_isControlled = true;
	}

	return m_dsv;
}

double Bus::GetDsv(void) {
	return m_dsv;
}

Bus::Violation Bus::GetStatus(void) {
	return m_status;
}

bool Bus::IsControlled(void) {
	return m_isControlled;
}

double Bus::GetACalc() {
	return m_aCalc;
}

double Bus::GetVCalc() {
	return m_vCalc;
}

double Bus::GetPCalc() {
	return m_pgCalc;
}
double Bus::GetQCalc() {
	return m_qgCalc;
}

void Bus::SetACalc(double aCalc) {
	m_aCalc = aCalc;
}

void Bus::SetVCalc(double vCalc) {
	m_vCalc = vCalc;
}

void Bus::SetPCalc(double pCalc) {
	m_pgCalc = pCalc;
}

void Bus::SetQCalc(double qCalc) {
	m_qgCalc = qCalc;
}

void Bus::SetOrdG(int ordG) {
	m_ordG = ordG;
}

int Bus::GetOrdG(void) const {
	return m_ordG;
}

void Bus::AddControl(double value) {
	m_control.push_back(value);
}

vec Bus::GetControl(void) const {
	return m_control;
}

void Bus::ClearControl(void) {
	m_control.clear();
}

double Bus::CalcPG(void) {
	double pg = this->GetBus().m_pc;
	std::vector<Bus*> neighbors = this->GetNeighbors();
	std::vector<Branch*> branches = this->GetBranches();
	int size = (int) branches.size();
	for (int j = 0; j < size; j++) {
		DBranch_t dataBranch = branches.at(j)->GetBranch();
		Bus* busM = neighbors.at(j);

		double theta_km = this->GetACalc() - busM->GetACalc();
		//Pg = Pg + (gkm(km)*V(slack)^2 - V(slack)*V(m)*(gkm(km)*cos(akm)+bkm(km)*sin(akm)));
		pg += (dataBranch.m_g * pow(this->GetVCalc(), 2) - this->GetVCalc() * busM->GetVCalc()
			  * (dataBranch.m_g * cos(theta_km) + dataBranch.m_b * sin(theta_km)));
		this->SetPCalc(pg);
	}
	return pg;

}

void Bus::SetVPreBusca(double v) {
	m_vPreBusca = v;
}

double Bus::GetVPreBusca(void) const {
	return m_vPreBusca;
}

void Bus::SetAPreBusca(double a) {
	m_aPreBusca = a;
}

double Bus::GetAPreBusca(void) const {
	return m_aPreBusca;
}

}
