/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef RUNTIME_EXCEPTION_H_
#define RUNTIME_EXCEPTION_H_

#include <exception>
#include <string>

using namespace std;

namespace exceptions {

class RuntimeException: public exception {
public:
	RuntimeException(string msg);
	virtual ~RuntimeException() throw();
	virtual const char* what() const throw ();
private:
	string m_msg;
};

} /* namespace exceptions */

#endif /* RUNTIME_EXCEPTION_H_ */
