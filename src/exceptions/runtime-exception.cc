/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#include "runtime-exception.h"

namespace exceptions {

RuntimeException::RuntimeException(string msg) {
	m_msg = msg;
}

RuntimeException::~RuntimeException() throw () {
}

const char* RuntimeException::what() const throw () {
	return m_msg.c_str();
}

} /* namespace exceptions */
