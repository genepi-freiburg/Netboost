#ifndef __functiontracker_hpp
#define __functiontracker_hpp

#include <iostream>
#include <string>
#include <sys/types.h>
#include <unistd.h>

//lonshy fix - using std 
using std::string;
using std::ostream;
using std::endl;
//end

struct functiontracker {
#ifdef FUNCTIONTRACKINGON
  string fn;
  ostream& out;

  functiontracker(string fn, ostream & out = cerr) :
    fn(fn),
    out(out) {
    out << "in " << fn << "(pid=" << getpid() << ")" << endl;
  }

  ~functiontracker() {
    out << "out " << fn << "(pid=" << getpid() << ")" << endl;
  }
#else
  functiontracker(string) {
  }

  functiontracker(string, ostream &) {
  }
#endif
};

#endif
