#ifndef __Readable_hpp
#define __Readable_hpp

#include <iostream>

template <typename tConcreteClass>
class tReadable {
};

template <typename tConcreteClass>
istream& operator>>(istream & is, tReadable<tConcreteClass> & r) {
  return static_cast<tConcreteClass &>(r).Read(is);
}

#endif
