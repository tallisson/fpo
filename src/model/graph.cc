/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */


#include "graph.h"

namespace model {

Graph::Graph() {
	m_numBus = 0;
	m_numBranch = 0;
	m_numPQ = 0;
	m_numPV = 0;
	m_numSlack = 0;
	m_posSlack = 0;
	m_slack = 0;
	m_contBusG = 0;
}

Graph::~Graph() {
	m_buses.clear();
	m_branches.clear();
	m_numBus = 0;
	m_numBranch = 0;
	m_numPQ = 0;
	m_numPV = 0;
	m_numSlack = 0;
}

void Graph::AddBus(Bus* bus) {
	int id = bus->GetBus().m_nin;
	m_buses.insert(std::pair<int, Bus*>(id, bus));
	m_numBus++;

	if (bus->GetBus().m_tipo == Bus::SLACK) {
		bus->SetType(Bus::SLACK);
		bus->SetVCalc(bus->GetBus().m_v);
		bus->SetACalc(bus->GetBus().m_ang);
		bus->SetPCalc(bus->GetBus().m_pg);
		bus->SetQCalc(bus->GetBus().m_qg);
		m_numSlack++;
		m_slack = bus;
		bus->SetOrdG(m_contBusG++);
	}

	if (bus->GetBus().m_tipo == Bus::GENERATION) {
		bus->SetType(Bus::GENERATION);
		bus->SetVCalc(bus->GetBus().m_v);
		bus->SetPCalc(bus->GetBus().m_pg);
		bus->SetQCalc(bus->GetBus().m_qg);
		m_numPV++;
		bus->SetOrdG(m_contBusG++);
	}

	if (bus->GetBus().m_tipo == Bus::LOAD
			|| bus->GetBus().m_tipo == Bus::LOSS_CONTROL_REACT) {
		if (bus->GetBus().m_tipo == 0) {
			bus->SetType(Bus::LOAD);
		} else if (bus->GetBus().m_tipo == 4) {
			bus->SetType(Bus::LOSS_CONTROL_REACT);
		}
		bus->SetPCalc(bus->GetBus().m_pg);
		bus->SetQCalc(bus->GetBus().m_qg);
		m_numPQ++;
	}

}

void Graph::Assoc(Branch* branch) {
	int idK = branch->GetBranch().m_nf;
	int idM = branch->GetBranch().m_ni;

	std::map<int, Bus*>::iterator itK, itM;
	Bus* busK = m_buses.find(idK)->second;
	Bus* busM = m_buses.find(idM)->second;

	busK->AddBranch(branch, busM);
	busM->AddBranch(branch, busK);

	m_branches.push_back(branch);
	m_numBranch++;
}

std::map<int, Bus*> Graph::GetBuses() const {
	return m_buses;
}

int Graph::GetNumBus() const {
	return m_numBus;
}

int Graph::GetNumBranch() const {
	return m_numBranch;
}

Bus*
Graph::GetBus(int idBus) {
	std::map<int, Bus*>::iterator it = m_buses.find(idBus);
	Bus* bus = 0;
	if (it != m_buses.end()) {
		bus = it->second;
	}

	return bus;
}

std::vector<Branch*> Graph::GetBranches() const {
	return m_branches;
}

int Graph::GetPosSlack(void) const {
	return m_posSlack;
}

void Graph::SetPosSlack(int posSlack) {
	m_posSlack = posSlack;
}

int Graph::GetNumPQ() const {
	return m_numPQ;
}

void Graph::SetNumPQ(int numPQ) {
	m_numPQ = numPQ;
}

int Graph::GetNumPV() const {
	return m_numPV;
}

void Graph::SetNumPV(int numPV) {
	m_numPV = numPV;
}

Bus* Graph::GetSlack(void) const {
	return m_slack;
}

void Graph::ClearContG(void) {
	m_contBusG = 0;
}

}
