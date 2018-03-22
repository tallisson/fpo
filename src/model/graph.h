/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "branch.h"
#include "bus.h"

#include <string>
#include <map>

namespace model {

class Graph {
private:
	int m_numBus;
	int m_numBranch;

	int m_numPQ;
	int m_numPV;
	int m_numSlack;
	int m_posSlack;
	int m_contBusG;

	std::map<int, Bus*> m_buses;
	std::vector<Branch*> m_branches;

	Bus* m_slack;
public:
	Graph();
	virtual ~Graph();

	void AddBus(Bus* bus);
	void Assoc(Branch* branch);

	int GetPosSlack(void) const;
	void SetPosSlack(int posSlack);

	std::map<int, Bus*> GetBuses() const;
	Bus* GetBus(int idBus);

	std::vector<Branch*> GetBranches() const;

	int GetNumBus() const;
	int GetNumBranch() const;

	int GetNumPQ() const;
	void SetNumPQ(int numPQ);

	void SetNumPV(int numPV);
	int GetNumPV() const;

	Bus* GetSlack(void) const;

	void ClearContG(void);
};

}

#endif  // GRAPH_H_
