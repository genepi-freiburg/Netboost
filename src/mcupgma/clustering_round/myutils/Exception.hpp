#ifndef __Exception_hpp
#define __Exception_hpp

#include <myutils/Writable.hpp>
#include <myutils/RCPtr.hpp>
#include <string>
#include <sstream>

class tException : public tWritable<tException> {
public:
  tException();

  tException(tException const &);

  ostream & Write(ostream& os) const;

  template<typename T>
  tException& operator<<(T const &);

protected:
  tException(char const * const type);
  
private:
  RCPtr<string> type;
  RCPtr<std::ostringstream> s;
};

inline
tException::tException() :
  type(new string("tException")),
  s(new std::ostringstream()) {
}

inline
tException::tException(tException const & other) :
  type(other.type),
  s(other.s) {
}

inline
tException::tException(char const * const type) :
  type(new string(type)),
  s(new std::ostringstream()) {
}

inline
ostream& tException::Write(ostream& os) const {
  return os << *type << ' ' << s->str();
}

template<typename T>
tException& tException::operator<<(T const & t) {
  *s << t;
  return *this;
}

#endif
