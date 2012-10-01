//----------------------------------------------------------------------------
//  File:    Slurper.cc
//  Purpose: Slurp tabular data from a text file.
//  Created: 03-May-2005 Harrison B. Prosper
//           20-Oct-2006 HBP Added mget method
//$Revision: 1.3 $ 
//----------------------------------------------------------------------------
#ifdef PROJECT_NAME
#include "PhysicsTools/LiteAnalysis/interface/Slurper.hpp"
#else
#include "Slurper.h"
#endif

#ifdef __WITH_CINT__
  ClassImp(Slurper)
#endif

#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
//----------------------------------------------------------------------------

using namespace std;

// Hide error within an anonymous namespace so that
// it is visible only within this programming unit
namespace {

  void error(string message)
  {
    cerr << "read ** error *** " << message << endl;
  }

  int wc(string filename)
  {
    int recordsize=1024;
    char record[1024];
    sprintf(record, "wc -l %s", filename.c_str());

    FILE* f = popen(record, "r");
    fread(record, 1, recordsize, f);
    pclose(f);

    istringstream inp(record);
    int nline;
    inp >> nline; 
    return nline;
  }
}

Slurper::Slurper()
  : _filename(""),
    _start(0),
    _count(0),
    _nrow(0),
    _size(0),
    _ok(true),
    _opened(false)
{}

Slurper::~Slurper() { close(); }

Slurper::Slurper(string filename, int start, int count)
  : _filename(filename),
    _start(start),
    _count(count),
    _nrow(0),
    _size(0),
    _ok(true),
    _opened(false)
{
  _stream = new ifstream(_filename.c_str());
  //  _stream.open(_filename.c_str());
  if ( ! _stream || ( _stream && !_stream->good()) ) 
    {
       string message = string("Slurper - unable to open ")+_filename;
       error(message);
      _ok = false;
    }
  else
    _opened = true;

  // Get number of entries in file

  _size = wc(_filename) - 1; 

  // Read header

  string line;
  getline(*_stream, line, '\n');
  istringstream inp(line);

  string name;
  _data.clear();
  _buffer.clear();
  _name.clear();
  int index=0;
  while ( inp >> name ) 
    {
      _data.push_back(0);
      _var[name] = index;
      _name.push_back(name);
      index++;
    }

  // Skip the first "start" lines

  for(int i=0; i < _start; i++)
    if ( !getline(*_stream, line, '\n') ) break;
}

bool
Slurper::read()
{
  _ok = true;

  // Read "count" rows if count > 0, otherwise read all lines

  string line;
  if ( !getline(*_stream, line, '\n') )
    {
      _ok = false;
      return _ok;
    }

  istringstream inp(line);
  for(int i=0; i < (int)_var.size(); i++) inp >> _data[i];

  _nrow++;

  if ( _count <= 0 )
    return true;
  else if (_nrow < _count )
    return true;
  else  
    return false;
}

int 
Slurper::entry() { return _nrow + _start; }

bool
Slurper::good() { return _ok; }

void
Slurper::close() 
{ 
  if (_opened)
    { 
      _stream->close();
      delete _stream;
      _stream = 0;
      _opened = false; 
    }
}

void
Slurper::rewind()
{
  _nrow = 0;
  _stream->seekg(ios_base::beg);
}

int
Slurper::size() { return _size; }

int
Slurper::entries() { return size(); }

bool
Slurper::present(string name) 
{ 
  return _var.find(name) != _var.end();
}

double
Slurper::get(string name) 
{
  if ( present(name) )
    return _data[_var[name]];
  else
    return 0;
}


std::vector<string>
Slurper::nget() 
{
    return _name;
}



std::vector<double>
Slurper::vget(std::vector<std::string>& names)
{
  if ( _buffer.size() < names.size() )
    {
      int nitem = names.size() - _buffer.size();
      for(int i=0; i < nitem; i++)
        _buffer.push_back(0);
    }
  for(int i=0; i < (int)names.size(); i++)
    {
      if ( present(names[i]) )
        _buffer[i] = _data[_var[names[i]]];
      else
        cout << "Slurper::vget -  name " << names[i] << " not found" << endl;
    }
  return _buffer;
}


void
Slurper::mget(std::map<std::string, double>& data)
{
  std::map<std::string, int>::iterator iter = _var.begin();
  while ( iter != _var.end() )
    {
      data[iter->first] = _data[_var[iter->first]];
      iter++;
    }
}
