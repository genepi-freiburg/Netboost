#ifndef __StreamFromFileNameGenerator_hpp
#define __StreamFromFileNameGenerator_hpp

#include <myutils/RCPtr.hpp>
#include <iostream>
#include <fstream>
#include <string>

inline
RCPtr<istream> InputStream(string const & filename) {
  if (filename == "-")
    return RCPtr<istream>(&cin, false);
//  else 
//    return new ifstream(filename.c_str());
// lonshy bug fix for silencing bad file
   else {
	ifstream * i = new ifstream(filename.c_str());
	if (!i->is_open()){
		delete i;
		i = NULL;
		throw tException() << "can't open input file:" << filename ;
	}
	return i;
   }


}

inline
RCPtr<ostream> OutputStream(string const & filename) {
  if (filename == "-")
    return RCPtr<ostream>(&cout, false);
//  else
//    return new ofstream(filename.c_str());
// lonshy bug fix for silencing bad file
    else {
	ofstream * o = new ofstream(filename.c_str());
	if (!o->is_open()){
		delete o;
		o = NULL;
		throw tException() << "can't open output file:" << filename ;
	}
	return o;
    }
}


#endif

