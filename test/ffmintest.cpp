/**********************************************************************
ffbfgs.cpp - Test BFGS minimization

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

Some portions Copyright (C) 2012 David C. Lonie
 
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

#include <string>
#include <cstdio>

#include <openbabel/math/vector3.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/obutil.h>

using namespace std;
using namespace OpenBabel;

#ifndef TESTDATADIR
#define TESTDATADIR="files/";
#endif

// Creates a "random" vector3 based on the value of seed. This is to ensure
// reproducibility.
const double DISPLACEMENT = 1.0;
inline vector3 generatePseudoRandomVector3(double &seed) {
  seed = fmod(seed, 1.0);
  double x, y, z;
  seed *= -2.4;
  seed += -0.8;
  x = seed;
  seed *= -0.5;
  seed += -1.8;
  y = seed;
  seed *= 1.5;
  seed += 2.3;
  z = seed;
  vector3 vec (x,y,z);
  vec.normalize();
  vec *= DISPLACEMENT;
//  cout << vec.x() << " " << vec.y() << " " << vec.z() << endl;
  return vec;
}

int main(int argc,char *argv[])
{
  // Define location of file formats for testing
#ifdef FORMATDIR
  char env[BUFF_SIZE];
  snprintf(env, BUFF_SIZE, "BABEL_LIBDIR=%s", FORMATDIR);
  putenv(env);
#endif

  string testdatadir = TESTDATADIR;
  string molecule_file = testdatadir + "mmff94-minima.mol2";

  OBMol refMol;
  OBConversion conv;
  conv.SetInAndOutFormats("mol2", "cml");
  if (!conv.ReadFile(&refMol, molecule_file)) {
    cerr << "Error reading file " << molecule_file.c_str() << endl;
    return 1;
  }

  // Displace coordinates
  OBMol randMol (refMol);
  double seed = 32.0;
  FOR_ATOMS_OF_MOL(a, randMol) {
    // Generate random vector
    vector3 disp = generatePseudoRandomVector3(seed);
    a->SetVector(a->GetVector() + disp);
  }

  // Optimize test molecule
  OBForceField *ff = OBForceField::FindForceField("MMFF94");
  if (ff == NULL) {
    cerr << "Cannot load force field!" << endl;
  }

  ff->SetLogFile(&cout);
  ff->SetLogLevel(OBFF_LOGLVL_NONE);
  ff->SetLineSearchType(LineSearchType::Newton2Num);
  ff->EnableCutOff(false);

  const int numSteps = 3000;
  const double econv = 1e-10;

  // Conjugate Gradients
  OBMol conjugateGradientsMol;
  ff->Setup(randMol);
  ff->ConjugateGradients(numSteps, econv);
  ff->GetCoordinates(conjugateGradientsMol);
  double conjugateGradientsEnergy = ff->Energy();

  // Steepest Descent
  OBMol steepestDescentMol;
  ff->Setup(randMol);
  ff->SteepestDescent(numSteps, econv);
  ff->GetCoordinates(steepestDescentMol);
  double steepestDescentEnergy = ff->Energy();

  // BFGS
  OBMol BFGSMol;
  ff->Setup(randMol);
  ff->SetLogLevel(OBFF_LOGLVL_LOW);
  ff->BFGS(numSteps, econv);
  ff->SetLogLevel(OBFF_LOGLVL_NONE);
  ff->GetCoordinates(BFGSMol);
  double BFGSEnergy = ff->Energy();

  conv.WriteFile(&refMol,                "/tmp/ffMinTest-ref.cml");
  conv.WriteFile(&randMol,               "/tmp/ffMinTest-rand.cml");
  conv.WriteFile(&steepestDescentMol,    "/tmp/ffMinTest-steepest.cml");
  conv.WriteFile(&conjugateGradientsMol, "/tmp/ffMinTest-conjugate.cml");
  conv.WriteFile(&BFGSMol,               "/tmp/ffMinTest-bfgs.cml");

  // Evaluate reference energy
  ff->Setup(refMol);
  double refEnergy = ff->Energy();

  // Compare energies:
  const double checkTol = 1.0;
  std::printf("%-20s %9.5f\n", "Reference:", refEnergy);
  std::printf("%-20s %9.5f %-7s\n",
              "SteepestDescent:",
              steepestDescentEnergy,
              fabs(steepestDescentEnergy - refEnergy) < checkTol
              ? "ok!" : "not ok!");
  std::printf("%-20s %9.5f %-7s\n",
              "ConjugateGradients:",
              conjugateGradientsEnergy,
              fabs(conjugateGradientsEnergy - refEnergy) < checkTol
              ? "ok!" : "not ok!");
  std::printf("%-20s %9.5f %-7s\n",
              "BFGS:",
              BFGSEnergy,
              fabs(BFGSEnergy - refEnergy) < checkTol
              ? "ok!" : "not ok!");

  return 0;
}

