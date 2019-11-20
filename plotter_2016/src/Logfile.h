////////////////////////
////  LOGFILE CLASS ////
////////////////////////

/*

Class is a container for the logfile
Essentially a container for the iostream
A little overkill, but this system will allow for more
complex logfiles to be made easily.  Also makes the 
code a little cleaner.  Who doesn't love clean code?

*/


#ifndef _LOGFILE_H_
#define _LOGFILE_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

using namespace std;

class Logfile {
 public:
  /// Constructor + Destructor
  //// if other is added, put copy constructor and operator=
  Logfile() { logfile.open("log.txt", ios::out);  }
  ~Logfile() { logfile.close(); }
  Logfile(string logname) { logfile.open(logname, ios::out); }


  /////Helper functions
  //// does the actual writing for the class
  void setHeader(vector<string>);
  void addLine(vector<string>);
  void setTrailer();

 private:
  ofstream logfile;
};

#endif
