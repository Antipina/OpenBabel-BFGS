/**********************************************************************
mol.cpp - Unit tests for Open Babel OBMol class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

string cmlfile = "cmltest/cs2a.cml";
string cmlfile_multi = "files/3d.head.2.cml";

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: mol" << endl;
      cout << " Unit tests for OBMol " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBMol \n";

  cout << "ok 1\n"; // for loading tests

  OBConversion obconv_first;
  if(!obconv_first.SetInFormat("CML"))
  {
    cout << "Bail out! Fail format isn't loaded!" << endl;
    return -1;
  }

  // Test using ReadFile to read from CML
  OpenBabel::OBMol obmol_first;
  if (obconv_first.ReadFile(&obmol_first, cmlfile))
          cout << "ok 2!" << endl;
  else
          cout << "not ok 2" << endl;


  
  OBConversion obconv;
  if(!obconv.SetInFormat("CML"))
  {
    cout << "Bail out! Fail format isn't loaded!" << endl;
    return -1;
  }
  // Test using ReadFile to read from multimol CML
  OpenBabel::OBMol obmol;
  if (obconv.ReadFile(&obmol, cmlfile_multi))
          cout << "ok 3!" << endl;
  else
          cout << "not ok 3" << endl;

  // Test reading the second and final molecule using Read
  if (obconv.Read(&obmol))
          cout << "ok 4!" << endl;
  else
          cout << "not ok 4" << endl;
  if (obmol.NumAtoms()==28)
          cout << "ok 5!" << endl;
  else
          cout << "not ok 5" << endl;
 
  // Test that there are no remaining molecules
  // (this test fails on Linux)
   if (!obconv.Read(&obmol))
          cout << "ok 6!" << endl;
   else
          cout << "not ok 6" << endl;

  // the total number of tests for "prove"
  // update when you add more tests!
  cout << "1..5\n";

  return 0;
}
