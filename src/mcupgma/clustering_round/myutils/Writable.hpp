#ifndef __Writable_hpp
#define __Writable_hpp

#include <iostream>

using namespace std;

template <typename tConcreteClass>
class tWritable {
};

template <typename tConcreteClass>
ostream& operator<<(ostream & os, tWritable<tConcreteClass> const & w) {
  return static_cast<tConcreteClass const &>(w).Write(os);
}

#endif
