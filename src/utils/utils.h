/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "../model/bus.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace model;

namespace utils {

struct sts
{
public:
  bool m_hasTape;
  double m_baseMVA;

  int m_nb;
  int m_npv;
  int m_npq;
  int m_nc;
  int m_ntap;
  int m_posSlack;

  std::vector<DBus_t> buses;
  std::vector<DBranch_t> branches;

  int BusFromNex (int nex);
};
typedef sts Sts_t;

class Utils
{
private:
	Sts_t m_s;
public:
	Utils();
	virtual ~Utils();

  Sts_t Read(std::string filename);

  static const double STEP = 0.00625;
  static const double TAPMIN_MIN = 0.75;
  static const double TAPMAX_MAX = 1.25;

  static const double TAPMIN = 0.90;
  static const double TAPMAX = 1.10;

  static const double VMIN_MIN = 0.65;
  static const double VMAX_MAX = 1.35;

  static const double VMIN = 0.90;
  static const double VMAX = 1.10;

  Sts_t GetSts(void) const;

  std::string Format (double v, int p, int maxS);
};

}
#endif /* UTILS_H_ */
