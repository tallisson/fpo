/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef LOAD_FLOW_H_
#define LOAD_FLOW_H_

#include "../model/graph.h"
#include "jacobian.h"
#include "mismatches.h"
#include "../utils/utils.h"
#include "../control/control.h"
#include "../control/q-control.h"
#include "../control/max-q.h"
#include "../control/vsf.h"

#include <string>
#include <sstream>

using namespace model;
using namespace utils;
using namespace control;


namespace load {

class LoadFlow {
private:
	Graph* m_graph;
	Control* m_qControl;
	Mismatch* m_mismatches;
	Jacobian* m_jac;
	VControl* m_vControl;

	Sts_t m_sts;

	double m_precision;
	double m_totalL;
	uint32_t m_iter;
	uint32_t m_maxIter;

	vec m_losses;
	vec m_b;

	void InitX0(void);

	void InitJ(void);

	void CalcLosses(void);

	bool m_verbose;
	bool m_hasQControl;
	std::string m_dir;
	std::string m_file;

	std::ostringstream os;

	int m_conv;
public:
	LoadFlow(void);
	virtual ~LoadFlow(void);

	Graph* GetGraph(void) const;

	Control* GetQControl(void) const;
	void SetQControl(Control* qControl);

	void SetVControl(VControl* vControl);

	void Prepare(std::string cdf);

	int Execute();

	void SetPrecision(double precision);

	void SetDir(std::string dir);

	void SetFile(std::string file);

	void SetQControl(bool hasQControl);

	void Reset(void);

	void SetVerbose(bool verbose);
	bool GetVerbose(void) const;

	mat ExecuteJ(void);

	int GetConv(void) const;

	Jacobian* GetJac(void) const;
};

}

#endif /* LOAD_FLOW */
