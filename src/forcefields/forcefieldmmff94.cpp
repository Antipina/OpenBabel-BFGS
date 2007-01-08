/*********************************************************************
forcefieldmmff94.cpp - Merck Molecular Force Field.
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <iomanip>
#include "forcefieldmmff94.h"

using namespace std;

namespace OpenBabel
{
  double OBFFBondCalculationMMFF94::Result()
  {
    double delta2;
    vector3 vab;
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    delta = rab - r0;
    delta2 = delta * delta;
 
    e = 143.9325f * 0.5f * kb * delta2 * (1.0f - 2.0f * delta + 7/12 * 4.0f * delta2);

    return e;
  }
  
  double OBForceFieldMMFF94::E_Bond()
  {
    vector<OBFFBondCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[100];

    energy = 0.0f;

    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "B O N D   S T R E T C H I N G" << endl << endl;
      *logos << "ATOM TYPES   FF    BOND       IDEAL       FORCE" << endl;
      *logos << " I    J     CLASS  LENGTH     LENGTH     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "------------------------------------------------------------------------" << endl;
    }
    
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); i++) {

      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).e);
        *logos  << logbuf << std::endl;
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM
      *logos << endl << "     TOTAL BOND STRETCHING ENERGY = " << energy << endl << endl;
    return energy;
  }
 
  double OBFFAngleCalculationMMFF94::Result()
  {
    double delta2;

    theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    
    delta = theta - theta0;
    delta2 = delta * delta;
    
    e = 0.043844f * 0.5f * ka * delta2 * (1.0f - 0.007f * delta);
     
    return e;
  }
  
  double OBForceFieldMMFF94::E_Angle()
  {
    vector<OBFFAngleCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[100];
 
    energy = 0.0f;

    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "A N G L E   B E N D I N G" << endl << endl;
      *logos << "ATOM TYPES        FF    VALENCE     IDEAL      FORCE" << endl;
      *logos << " I    J    K     CLASS   ANGLE      ANGLE     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "-----------------------------------------------------------------------------" << endl;
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); i++) {

      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d      %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).at, (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).e);
        *logos  << logbuf << endl;
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM
      *logos << endl << "     TOTAL ANGLE BENDING ENERGY = " << energy << endl << endl;
    return energy;
  }
  
  double OBFFStrBndCalculationMMFF94::Result()
  {
    vector3 vab, vbc;

    vab = a->GetVector() - b->GetVector();
    vbc = b->GetVector() - c->GetVector();
    
    theta = a->GetAngle(b->GetIdx(), c->GetIdx());
    rab = vab.length();
    rbc = vbc.length();
    
    delta_theta = theta - theta0;
    delta_rab = rab - rab0;
    delta_rbc = rbc - rbc0;

    e = 2.51210f * (kbaABC * delta_rab + kbaCBA * delta_rbc) * delta_theta;

    return e;
  }
  
  double OBForceFieldMMFF94::E_StrBnd() 
  {
    vector<OBFFStrBndCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[100];

    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos  << std::endl << "S T R E T C H   B E N D I N G" << std::endl << std::endl;
      *logos << "ATOM TYPES        FF    VALENCE     DELTA        FORCE CONSTANT " << std::endl;
      *logos << " I    J    K     CLASS   ANGLE      ANGLE        I J        J K      ENERGY" << std::endl;
      *logos << "---------------------------------------------------------------------------" << std::endl;
    }
    
    for (i = _strbndcalculations.begin(); i != _strbndcalculations.end(); i++) {

      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d     %2d   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), 
                (*i).sbt, (*i).theta, (*i).delta_theta, (*i).kbaABC, (*i).kbaCBA, (*i).e);
        *logos  << logbuf << endl;
      }
    }
	
    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL STRETCH BENDING ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
 
  int OBForceFieldMMFF94::GetElementRow(OBAtom *atom)
  {
    int row;
    
    row = 0;

    if (atom->GetAtomicNum() > 2)
      row++;
    if (atom->GetAtomicNum() > 10)
      row++;
    if (atom->GetAtomicNum() > 18)
      row++;
    if (atom->GetAtomicNum() > 36)
      row++;
    if (atom->GetAtomicNum() > 54)
      row++;
    if (atom->GetAtomicNum() > 86)
      row++;
    
    return row;
  }

  double OBFFTorsionCalculationMMFF94::Result()
  {
    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;

    tor = CalcTorsionAngle(a->GetVector(), b->GetVector(), c->GetVector(), d->GetVector());

    cosine = cos(DEG_TO_RAD * tor);
    cosine2 = cos(DEG_TO_RAD * 2 * tor);
    cosine3 = cos(DEG_TO_RAD * 3 * tor);
      
    phi1 = 1.0f + cosine;
    phi2 = 1.0f - cosine2;
    phi3 = 1.0f + cosine3;
    
    e = 0.5f * (v1 * phi1 + v2 * phi2 + v3 * phi3);
 
    return e;
  }
  
  double OBForceFieldMMFF94::E_Torsion() 
  {
    vector<OBFFTorsionCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[150];
    
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "T O R S I O N A L" << std::endl << std::endl;
      *logos << "ATOM TYPES             FF     TORSION       FORCE CONSTANT" << std::endl;
      *logos << " I    J    K    L     CLASS    ANGLE         V1   V2   V3     ENERGY" << std::endl;
      *logos << "--------------------------------------------------------------------" << std::endl;
    }

    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); i++) {

      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      %d   %8.3f   %6.3f   %6.3f   %6.3f   %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), atoi((*i).d->GetType()), 
                (*i).tt, (*i).tor, (*i).v1, (*i).v2, (*i).v3, (*i).e);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL TORSIONAL ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
 
  //
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /
  //  c
  //
  double OBFFOOPCalculationMMFF94::Result()
  {
    double angle2;

    angle = Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
    angle2 = angle * angle;
    
    e = 0.043844f * 0.5f * koop * angle2;

    return e;
  }

  double OBForceFieldMMFF94::E_OOP() 
  {
    vector<OBFFOOPCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[100];

    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "O U T - O F - P L A N E   B E N D I N G" << std::endl << std::endl;
      *logos << "ATOM TYPES             FF       OOP     FORCE " << std::endl;
      *logos << " I    J    K    L     CLASS    ANGLE   CONSTANT     ENERGY" << std::endl;
      *logos << "----------------------------------------------------------" << std::endl;
    }

    for (i = _oopcalculations.begin(); i != _oopcalculations.end(); i++) {
      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), atoi((*i).c->GetType()), atoi((*i).d->GetType()), 
                (*i).angle, (*i).koop, (*i).e);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_HIGH
     *logos << std::endl << "     TOTAL OUT-OF-PLANE BENDING ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }
 
  double OBFFVDWCalculationMMFF94::Result()
  {
      double rab2, rab3, rab4, rab5, rab6, rab7;
      double erep2, erep3, erep4, erep5, erep6, erep7;

      rab = a->GetDistance(b);
      rab2 = rab * rab;
      rab4 = rab2 * rab2;
      rab6 = rab4 * rab2;
      rab7 = rab6 * rab;

      erep = (1.07f * R_AB) / (rab + 0.07f * R_AB); //***
      erep2 = erep * erep;
      erep4 = erep2 * erep2;
      erep6 = erep4 * erep2;
      erep7 = erep6 * erep;
      
      eattr = (((1.12f * R_AB7) / (rab7 + 0.12f * R_AB7)) - 2.0f);
      
      e = escale * epsilon * erep7 * eattr;
      
      return e;
  }

  double OBForceFieldMMFF94::E_VDW()
  {
    vector<OBFFVDWCalculationMMFF94>::iterator i;
    double energy;
    char logbuf[100];

    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "V A N   D E R   W A A L S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J        Rij       R*IJ    EPSILON    E_REP     E_ATTR    ENERGY" << std::endl;
      *logos << "----------------------------------------------------------------------" << std::endl;
      //          XX   XX     -000.000  -000.000  -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); i++) {
      
      energy += (*i).Result();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%2d   %2d     %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f", atoi((*i).a->GetType()), atoi((*i).b->GetType()), 
	        (*i).rab, (*i).R_AB, (*i).epsilon, (*i).erep, (*i).eattr, (*i).e);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL VAN DER WAALS ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }

  double OBForceFieldMMFF94::E_Electrostatic()
  {
    OBAtom *a, *b, *c, *d;
    double e, energy;
    char logbuf[100];
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "E L E C T R O S T A T I C   I N T E R A C T I O N S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J        Rij       R*IJ    EPSILON    E_REP     E_ATTR    ENERGY" << std::endl;
      *logos << "----------------------------------------------------------------------" << std::endl;
      //            XX   XX     -000.000  -000.000  -000.000  -000.000  -000.000  -000.000
    }

    energy = 0.0f;
    
    FOR_PAIRS_OF_MOL(p, _mol) {
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "     TOTAL VAN DER WAALS ENERGY = " << energy << std::endl << std::endl;
    return energy;
  }

  //
  // OBForceFieldMM2 member functions
  //
  //***********************************************
  //Make a global instance
  OBForceFieldMMFF94 theForceFieldMMFF94("MMFF94",true);
  //***********************************************

  OBForceFieldMMFF94::~OBForceFieldMMFF94()
  {
  }

  OBForceFieldMMFF94 &OBForceFieldMMFF94::operator=(OBForceFieldMMFF94 &src)
  {
    _mol = src._mol;
  }

  bool OBForceFieldMMFF94::Setup(OBMol &mol)
  {
    UnsetEnergyCalculated();

    _mol = mol;
    SetMMFFTypes();
    SetupCalculations();
    //CalcCharges();
  }
 
  bool OBForceFieldMMFF94::ParseParamFile()
  {
    ParseParamProp();
    ParseParamBond();
    ParseParamBndk();
    ParseParamAngle();
    ParseParamStrBnd();
    ParseParamDfsb();
    ParseParamOOP();
    ParseParamTorsion();
    ParseParamVDW();
    ParseParamCharge();
  }
  
  bool OBForceFieldMMFF94::ParseParamBond()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffbond.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffbond.par";
    buffer2 += "mmffbond.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.ipar5 = atoi(vs[0].c_str());  // FF class
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.dpar1 = atof(vs[3].c_str());  // kb
	parameter.dpar2 = atof(vs[4].c_str());  // r0
	_ffbondparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamBndk()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffbndk.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffbndk.par";
    buffer2 += "mmffbndk.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.dpar1 = atof(vs[3].c_str());  // r0-ref
	parameter.dpar2 = atof(vs[4].c_str());  // kb-ref
	_ffbndkparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamAngle()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffang.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffang.par";
    buffer2 += "mmffang.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.ipar5 = atoi(vs[0].c_str());  // FF class
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.c = atoi(vs[3].c_str());
	parameter.dpar1 = atof(vs[4].c_str());  // ka
	parameter.dpar2 = atof(vs[5].c_str());  // theta0
	_ffangleparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
   
  bool OBForceFieldMMFF94::ParseParamStrBnd()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffstbn.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffstbn.par";
    buffer2 += "mmffstbn.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.ipar5 = atoi(vs[0].c_str());  // FF class
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.c = atoi(vs[3].c_str());
	parameter.dpar1 = atof(vs[4].c_str());  // kbaIJK
	parameter.dpar2 = atof(vs[5].c_str());  // kbaKJI
	_ffstrbndparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamDfsb()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffdfsb.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffdfsb.par";
    buffer2 += "mmffdfsb.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[0].c_str());
	parameter.b = atoi(vs[1].c_str());
	parameter.c = atoi(vs[2].c_str());
	parameter.dpar1 = atof(vs[3].c_str());  // kbaIJK
	parameter.dpar2 = atof(vs[4].c_str());  // kbaKJI
	_ffdfsbparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamOOP()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffoop.par";
    buffer2 += "mmffoop.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[0].c_str());
	parameter.b = atoi(vs[1].c_str());
	parameter.c = atoi(vs[2].c_str());
	parameter.d = atoi(vs[3].c_str());
	parameter.dpar1 = atof(vs[4].c_str());  // koop
	_ffoopparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldMMFF94::ParseParamTorsion()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmfftor.par";
    buffer2 += "mmfftor.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.ipar5 = atoi(vs[0].c_str());  // FF class
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.c = atoi(vs[3].c_str());
	parameter.d = atoi(vs[4].c_str());
	parameter.dpar1 = atof(vs[5].c_str());  // v1
	parameter.dpar2 = atof(vs[6].c_str());  // v2
	parameter.dpar3 = atof(vs[7].c_str());  // v3
	_fftorsionparams.push_back(parameter);
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
   
  bool OBForceFieldMMFF94::ParseParamVDW()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffvdw.par";
    buffer2 += "mmffvdw.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[0].c_str());
	parameter.dpar1 = atof(vs[1].c_str());  // alpha-i
	parameter.dpar2 = atof(vs[2].c_str());  // N-i
	parameter.dpar3 = atof(vs[3].c_str());  // A-i
	parameter.dpar4 = atof(vs[4].c_str());  // G-i
        if (EQn(vs[5].c_str(), "-", 1))
	  parameter.ipar1 = 0;
        else if (EQn(vs[5].c_str(), "D", 1))
	  parameter.ipar1 = 1;  // hydrogen bond donor
        else if (EQn(vs[5].c_str(), "A", 1))
	  parameter.ipar1 = 2;  // hydrogen bond acceptor
	_ffvdwparams.push_back(parameter);
    }

    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
   
  bool OBForceFieldMMFF94::ParseParamCharge()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffchg.par";
    buffer2 += "mmffchg.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[1].c_str());
	parameter.b = atoi(vs[2].c_str());
	parameter.dpar1 = atof(vs[3].c_str());  // bci
	_ffchgparams.push_back(parameter);
    }

    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
 
  bool OBForceFieldMMFF94::ParseParamPbci()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;
    
    // open data/mmffoop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffpbci.par";
    buffer2 += "mmffpbci.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	parameter.clear();
	parameter.a = atoi(vs[1].c_str());
	parameter.dpar1 = atof(vs[2].c_str());  // pbci
	parameter.dpar2 = atof(vs[3].c_str());  // fcadj
	_ffpbciparams.push_back(parameter);
    }

    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
 
  bool OBForceFieldMMFF94::ParseParamProp()
  {
    vector<string> vs;
    char buffer[80];
    
    // open data/mmffprop.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffprop.par";
    buffer2 += "mmffprop.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
	tokenize(vs, buffer);

	if (atoi(vs[6].c_str()) == 1)
          _arom.push_back(atoi(vs[0].c_str()));

	if (atoi(vs[7].c_str()) == 1)
          _lin.push_back(atoi(vs[0].c_str()));

	if (atoi(vs[8].c_str()) == 1)
          _sbmb.push_back(atoi(vs[0].c_str()));
    }

    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }

  bool OBForceFieldMMFF94::SetMMFFTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/mmffsymb.par
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "mmffsymb.par";
    buffer2 += "mmffsymb.par";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
    {
      ifs2.open(buffer2.c_str());
      ifsP = &ifs2;
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << std::endl << "A T O M   T Y P E S" << std::endl << std::endl;
    
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[0]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[1]));
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
          return false;
        }
 
        for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
          if (i->first->Match(_mol)) {
            _mlist = i->first->GetMapList();
            for (j = _mlist.begin();j != _mlist.end();++j) {
              _mol.GetAtom((*j)[0])->SetType(i->second);
	    }
          }
        }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      *logos << "IDX\tTYPE" << std::endl;
      FOR_ATOMS_OF_MOL (a, _mol)
        *logos << a->GetIdx() << "\t" << a->GetType() << std::endl;
    }
  }
  
  bool OBForceFieldMMFF94::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    bool found;

    // 
    // Bond Calculations
    //
    OBFFBondCalculationMMFF94 bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = GetBondType(a, b);
      
      parameter = GetParameterMMFF94(bondtype, atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
      if (parameter == NULL) {
        
        //parameter = GetParameter(atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbndkparams);
	if (parameter == NULL) {
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters", obError);
          exit(1);
	}
      } else {
        bondcalc.a = a;
	bondcalc.b = b;
        bondcalc.kb = parameter->dpar1;
        bondcalc.r0 = parameter->dpar2;
        bondcalc.bt = bondtype;

	_bondcalculations.push_back(bondcalc);
      }
    }

    //
    // Angle Calculations
    //
    OBFFAngleCalculationMMFF94 anglecalc;
    int angletype;
 
    _anglecalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      angletype = GetAngleType(a, b, c);

      parameter = GetParameterMMFF94(angletype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffangleparams);
      if (parameter == NULL) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all angle parameters", obError);
        exit(1);
      }
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;
      anglecalc.ka   = parameter->dpar1;
      anglecalc.theta0 = parameter->dpar2;
      anglecalc.at = angletype;
      
      _anglecalculations.push_back(anglecalc);
    }

    //
    // StrBnd Calculations
    //
    OBFFStrBndCalculationMMFF94 strbndcalc;
    int strbndtype, bondtype1, bondtype2;

    _strbndcalculations.clear();

    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);

      if (HasLinSet(atoi(b->GetType())))
        continue;

      strbndtype = GetStrBndType(a, b, c);
      bondtype1 = GetBondType(a, b);
      bondtype2 = GetBondType(b, c);
      angletype = GetAngleType(a, b, c);
      
      parameter = GetParameterMMFF94(strbndtype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffstrbndparams);
      if (parameter == NULL) {
        int rowa, rowb, rowc;

        rowa = GetElementRow(a);
        rowb = GetElementRow(b);
        rowc = GetElementRow(c);
        
	parameter = GetParameter(rowa, rowb, rowc, 0, _ffdfsbparams);
        if (parameter == NULL) {
          obErrorLog.ThrowError(__FUNCTION__, "Could not find all stretch-bend parameters", obError);
          exit(1);
	}

        if (rowa == parameter->a) {
          strbndcalc.kbaABC = parameter->dpar1;
          strbndcalc.kbaCBA = parameter->dpar2;
        } else {
          strbndcalc.kbaABC = parameter->dpar2;
          strbndcalc.kbaCBA = parameter->dpar1;
        }
      } else {
        if (atoi(a->GetType()) == parameter->a) {
          strbndcalc.kbaABC = parameter->dpar1;
          strbndcalc.kbaCBA = parameter->dpar2;
        } else {
          strbndcalc.kbaABC = parameter->dpar2;
          strbndcalc.kbaCBA = parameter->dpar1;
        }
      }

      parameter = GetParameterMMFF94(angletype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), 0, _ffangleparams);
      strbndcalc.theta0 = parameter->dpar2;
      
      parameter = GetParameterMMFF94(bondtype1, atoi(a->GetType()), atoi(b->GetType()), 0, 0, _ffbondparams);
      strbndcalc.rab0 = parameter->dpar2;
      
      parameter = GetParameterMMFF94(bondtype2, atoi(b->GetType()), atoi(c->GetType()), 0, 0, _ffbondparams);
      strbndcalc.rbc0 = parameter->dpar2;

      strbndcalc.a = a;
      strbndcalc.b = b;
      strbndcalc.c = c;
      strbndcalc.sbt = strbndtype;

      _strbndcalculations.push_back(strbndcalc);
    }

    //
    // Torsion Calculations
    //
    OBFFTorsionCalculationMMFF94 torsioncalc;
    int torsiontype;

    _torsioncalculations.clear();
 
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      torsiontype = GetTorsionType(a, b, c, d);
      
      parameter = GetParameterMMFF94(torsiontype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);

      found = false;
      parameter = GetParameterMMFF94(torsiontype, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()), _fftorsionparams);
      if (parameter == NULL) {
        for (int idx=0; idx < _fftorsionparams.size(); idx++) {  // *-XX-XX-XX 
          if (((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) && (torsiontype == _fftorsionparams[idx].ipar5) &&
	       (atoi(c->GetType()) == _fftorsionparams[idx].c) && (atoi(d->GetType()) == _fftorsionparams[idx].d)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].c) && (torsiontype == _fftorsionparams[idx].ipar5) &&
	       (atoi(c->GetType()) == _fftorsionparams[idx].b) && (atoi(d->GetType()) == _fftorsionparams[idx].d)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) && (torsiontype == _fftorsionparams[idx].ipar5) &&
	       (atoi(c->GetType()) == _fftorsionparams[idx].c) && (atoi(d->GetType()) == _fftorsionparams[idx].a)) ||
              ((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].c) && (torsiontype == _fftorsionparams[idx].ipar5) &&
	       (atoi(c->GetType()) == _fftorsionparams[idx].b) && (atoi(d->GetType()) == _fftorsionparams[idx].a)))
	  {
            torsioncalc.v1 = _fftorsionparams[idx].dpar1;
            torsioncalc.v2 = _fftorsionparams[idx].dpar2;
            torsioncalc.v3 = _fftorsionparams[idx].dpar3;
	    found = true;  
	  }
        }

        if (!found)
	  for (int idx=0; idx < _fftorsionparams.size(); idx++) {  // *-XX-XX-*
            if (((_fftorsionparams[idx].a == 0) && (atoi(b->GetType()) == _fftorsionparams[idx].b) &&  (torsiontype == _fftorsionparams[idx].ipar5) &&
	         (atoi(c->GetType()) == _fftorsionparams[idx].c) && (_fftorsionparams[idx].d == 0)) ||
                ((_fftorsionparams[idx].a == 0) && (atoi(c->GetType()) == _fftorsionparams[idx].b) &&  (torsiontype == _fftorsionparams[idx].ipar5) &&
	         (atoi(b->GetType()) == _fftorsionparams[idx].c) && (_fftorsionparams[idx].d == 0)))
	    {
              torsioncalc.v1 = _fftorsionparams[idx].dpar1;
              torsioncalc.v2 = _fftorsionparams[idx].dpar2;
              torsioncalc.v3 = _fftorsionparams[idx].dpar3;
	      found = true;  
	    }
          }

	if (!found) {
	  obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
          exit(1);
	}
      } else {
        torsioncalc.v1 = parameter->dpar1;
        torsioncalc.v2 = parameter->dpar2;
        torsioncalc.v3 = parameter->dpar3;
      }
      
      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
      torsioncalc.tt = torsiontype;

      _torsioncalculations.push_back(torsioncalc);
    }

    //
    // Out-Of-Plane Calculations
    //
    OBFFOOPCalculationMMFF94 oopcalc;

    _oopcalculations.clear();
 
    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      found = false;

      for (int idx=0; idx < _ffoopparams.size(); idx++) {
        if (atoi(b->GetType()) == _ffoopparams[idx].b) {
          a = NULL;
          c = NULL;
          d = NULL;

	  FOR_NBORS_OF_ATOM(nbr, b) {
	    if (a ==NULL)
	      a = (OBAtom*) &*nbr;
	    else if (c == NULL)
	      c = (OBAtom*) &*nbr;
	    else
	      d = (OBAtom*) &*nbr;
	  }
	  

	  if (((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
	      ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(d->GetType()) == _ffoopparams[idx].d)) ||
	      ((atoi(c->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
	      ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(c->GetType()) == _ffoopparams[idx].c) && (atoi(a->GetType()) == _ffoopparams[idx].d)) ||
	      ((atoi(a->GetType()) == _ffoopparams[idx].a) && (atoi(d->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)) ||
	      ((atoi(d->GetType()) == _ffoopparams[idx].a) && (atoi(a->GetType()) == _ffoopparams[idx].c) && (atoi(c->GetType()) == _ffoopparams[idx].d)))
	  {
	    found = true;

	    oopcalc.koop = _ffoopparams[idx].dpar1;
            
	    // A-B-CD || C-B-AD  PLANE = ABC
            oopcalc.a = a;
            oopcalc.b = b;
            oopcalc.c = c;
            oopcalc.d = d;
	    
	    _oopcalculations.push_back(oopcalc);

	    // C-B-DA || D-B-CA  PLANE BCD
	    oopcalc.a = d;
            oopcalc.d = a;
	
	    _oopcalculations.push_back(oopcalc);
            
	    // A-B-DC || D-B-AC  PLANE ABD
	    oopcalc.a = a;
	    oopcalc.c = d;
            oopcalc.d = c;
	    
	    _oopcalculations.push_back(oopcalc);
	  }

	  if ((_ffoopparams[idx].a == 0) && (_ffoopparams[idx].c == 0) && (_ffoopparams[idx].d == 0) && !found) // *-XX-*-*
	  {
	    oopcalc.koop = _ffoopparams[idx].dpar1;
	    
	    // A-B-CD || C-B-AD  PLANE = ABC
            oopcalc.a = a;
            oopcalc.b = b;
            oopcalc.c = c;
            oopcalc.d = d;
	    
	    _oopcalculations.push_back(oopcalc);

	    // C-B-DA || D-B-CA  PLANE BCD
	    oopcalc.a = d;
            oopcalc.d = a;
	
	    _oopcalculations.push_back(oopcalc);
            
	    // A-B-DC || D-B-AC  PLANE ABD
	    oopcalc.a = a;
	    oopcalc.c = d;
            oopcalc.d = c;
	    
	    _oopcalculations.push_back(oopcalc);
	  }
        }
      }
    }

    // 
    // VDW Calculations
    //
    OBFFVDWCalculationMMFF94 vdwcalc;

    _vdwcalculations.clear();
 
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      OBFFParameter *parameter_a, *parameter_b;
      parameter_a = GetParameter(atoi(a->GetType()), 0, 0, 0, _ffvdwparams);
      parameter_b = GetParameter(atoi(b->GetType()), 0, 0, 0, _ffvdwparams);
      if ((parameter_a == NULL) || (parameter_b == NULL)) {
        obErrorLog.ThrowError(__FUNCTION__, "Could not find all Van der Waals parameters", obError);
        exit(1);
      }
      
      vdwcalc.a = a;
      vdwcalc.alpha_a = parameter_a->dpar1;
      vdwcalc.Na = parameter_a->dpar2;
      vdwcalc.Aa = parameter_a->dpar3;
      vdwcalc.Ga = parameter_a->dpar4;
      vdwcalc.aDA = parameter_a->ipar1;
      
      vdwcalc.b = b;
      vdwcalc.alpha_b = parameter_b->dpar1;
      vdwcalc.Nb = parameter_b->dpar2;
      vdwcalc.Ab = parameter_b->dpar3;
      vdwcalc.Gb = parameter_b->dpar4;
      vdwcalc.bDA = parameter_b->ipar1;
      
      //these calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      double R_AA, R_BB, R_AB6, g_AB, g_AB2;
      double R_AB2, R_AB4, R_AB7, sqrt_a, sqrt_b;
 
      R_AA = vdwcalc.Aa * pow(vdwcalc.alpha_a, 0.25f);
      R_BB = vdwcalc.Ab * pow(vdwcalc.alpha_b, 0.25f);
      
      vdwcalc.escale = 1.0f;
      if (vdwcalc.aDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5f * (R_AA + R_BB);
	if (vdwcalc.bDA == 2) { // hydrogen bond acceptor
	  vdwcalc.R_AB = 0.8f * vdwcalc.R_AB;
          vdwcalc.escale = 0.5f;
	}
      } else if (vdwcalc.bDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5f * (R_AA + R_BB);
	if (vdwcalc.aDA == 2) { // hydrogen bond acceptor
	  vdwcalc.R_AB = 0.8f * vdwcalc.R_AB;
          vdwcalc.escale = 0.5f;
	}
      } else {
        g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
        g_AB2 = g_AB * g_AB;
        vdwcalc.R_AB =  0.5f * (R_AA + R_BB) * (1.0f + 0.2f * (1.0f - exp(-12.0f * g_AB2)));
      }
      
      R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
      R_AB4 = R_AB2 * R_AB2;
      R_AB6 = R_AB4 * R_AB2;
      vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;

      sqrt_a = sqrt(vdwcalc.alpha_a / vdwcalc.Na);
      sqrt_b = sqrt(vdwcalc.alpha_b / vdwcalc.Nb);
      vdwcalc.epsilon = (181.16f * vdwcalc.Ga * vdwcalc.Gb * vdwcalc.alpha_a * vdwcalc.alpha_b) / (sqrt_a + sqrt_b) * (1.0f / R_AB6);
 
      _vdwcalculations.push_back(vdwcalc);
    }

  }

  bool OBForceFieldMMFF94::CalcCharges()
  {
    OBFFParameter *parameter_a, *parameter_b, *parameter_ab;
    OBAtom *a;
    double q0i;
    
    FOR_ATOMS_OF_MOL (atom, _mol) {
      a = (OBAtom*) &*atom;
      parameter_a = GetParameter(atoi(atom->GetType()), 0, 0, 0, _ffpbciparams);

        if (atom->IsCarboxylOxygen())
            q0i = -0.500;
        else if (atom->IsPhosphateOxygen() && atom->GetHvyValence() == 1)
            q0i = -0.666;
        else if (atom->IsSulfateOxygen())
            q0i = -0.500;
 
      *logos << "fcharge = " << q0i << std::endl;
      //FOR_NBORS_OF_ATOM (nbr, a) {
      // q =  
      //}
    }
  }

  double OBForceFieldMMFF94::Energy()
  {
    double energy;

    SetEnergyCalculated();
    
    energy = E_Bond();
    energy += E_Angle();
    energy += E_StrBnd();
    energy += E_OOP();
    energy += E_Torsion();
    energy += E_VDW();
    //energy += E_Electrostatic();

    return energy;
  }

  // used to validate the implementation
  bool OBForceFieldMMFF94::Validate ()
  {
    OBConversion conv;
    OBFormat *format_in = conv.FindFormat("mol2");
    vector<string> vs;
    vector<int> types;
    char buffer[150], logbuf[100];
    bool molfound, atomfound;
    double etot, ebond, eangle, eoop, estbn, etor, evdw, eeq;
    int n;

    if (!format_in || !conv.SetInFormat(format_in)) {
      obErrorLog.ThrowError(__FUNCTION__, "Could not set mol2 input format", obError);
      exit (-1);
    }

    ifstream ifs, ifs2;
    ofstream ofs;

    ifs.open("MMFF94_dative.mol2");
    if (!ifs) {
      obErrorLog.ThrowError(__FUNCTION__, "Could not open ./MMFF94_dative.mol2", obError);
      exit (-1);
    }
 
    ifs2.open("MMFF94_opti.log");
    if (!ifs2) {
      obErrorLog.ThrowError(__FUNCTION__, "Coulg not open ./MMFF_opti.log", obError);
      exit(1);
    }
    
    ofs.open("MMFF94_openbabel.log");
    if (!ofs) {
      obErrorLog.ThrowError(__FUNCTION__, "Coulg not open ./MMFF_openbabel.log", obError);
      exit(1);
    }
    
    SetLogFile(&ofs);
   
    for (unsigned int c=1;; c++) {
      _mol.Clear();
      types.clear();
      if (!conv.Read(&_mol, &ifs))
        break;
      if (_mol.Empty())
        break;
      
      UnsetEnergyCalculated();
      SetMMFFTypes();
      SetupCalculations();
      
      n = _mol.NumAtoms() / 4;
      molfound = false;
      atomfound = false;
      while (ifs2.getline(buffer, 150)) {
	tokenize(vs, buffer);
	if (vs.size() < 2)
	  continue;
	
        string str(buffer);
	if (string::npos != str.find(_mol.GetTitle(),0))
	  molfound = true;

	if (atomfound) {
          if (n) {
	    types.push_back(atoi(vs[2].c_str()));
	    types.push_back(atoi(vs[5].c_str()));
	    types.push_back(atoi(vs[8].c_str()));
	    types.push_back(atoi(vs[11].c_str()));
	  } else {
	    if (vs.size() > 2)
	      types.push_back(atoi(vs[2].c_str()));
	    if (vs.size() > 5)
	      types.push_back(atoi(vs[5].c_str()));
	    if (vs.size() > 8)
	      types.push_back(atoi(vs[8].c_str()));
	    
	    atomfound = false;
	  }
	  n--;
        }
        
	if (molfound && EQn(buffer, " ATOM NAME  TYPE", 16))
	  atomfound = true;
	
	if (molfound) {
	  if (EQn(buffer, " Total ENERGY", 13))
	    etot = atof(vs[3].c_str());
	  if (EQn(buffer, " Bond Stretching", 16))
	    ebond = atof(vs[2].c_str());
	  if (EQn(buffer, " Angle Bending", 14))
	    eangle = atof(vs[2].c_str());
	  if (EQn(buffer, " Out-of-Plane Bending", 21))
	    eoop = atof(vs[2].c_str());
	  if (EQn(buffer, " Stretch-Bend", 13))
	    estbn = atof(vs[1].c_str());
	  if (EQn(buffer, "     Total Torsion", 18))
	    etor = atof(vs[2].c_str());
	  if (EQn(buffer, "     Net vdW", 12))
	    evdw = atof(vs[2].c_str());
	  if (EQn(buffer, " Electrostatic", 14)) {
	    eeq = atof(vs[1].c_str());
	    break;
	  }
	}
      } // while (getline)
      
      vector<int>::iterator i;
      int ni;
    
      cout << "--------------------------------------------------------------------------------" << endl;
      cout << "                                                                                " << endl;
      cout << "  MOLECULE : " << _mol.GetTitle() << endl;
      cout << "                                                                                " << endl;
      cout << "IDX  OB_TYPE  LOG_TYPE      RESULT                                              " << endl;
      cout << "----------------------------------                                              " << endl;
 
      ni = 1;
      for (i = types.begin(); i != types.end();i++) {
        if (ni > _mol.NumAtoms())
          continue;

        if (atoi(_mol.GetAtom(ni)->GetType()) == (*i))
          sprintf(logbuf, "%2d    %3d      %3d          PASSED", _mol.GetAtom(ni)->GetIdx(), atoi(_mol.GetAtom(ni)->GetType()), *i);
        else
          sprintf(logbuf, "%2d    %3d      %3d      XXX FAILED XXX", _mol.GetAtom(ni)->GetIdx(), atoi(_mol.GetAtom(ni)->GetType()), *i);
      
        cout << logbuf << endl;
        
        ni++;
      }

      double delta, err;
      cout << endl;
      cout << "TERM                   OB ENERGY    LOG ENERGY      DELTA     REL.ERR (%)" << endl;
      cout << "-------------------------------------------------------------------------" << endl;
    
      delta = (E_Bond() - ebond);
      err = delta / ebond * 100;
      sprintf(logbuf, "Bond Stretching        %8.3f      %8.3f     %8.3f     %5.2f", E_Bond(), ebond, delta, err);
      cout << logbuf << endl;
    
      delta = (E_Angle() - eangle);
      err = delta / eangle * 100;
      sprintf(logbuf, "Angle Bending          %8.3f      %8.3f     %8.3f     %5.2f", E_Angle(), eangle, delta, err);
      cout << logbuf << endl;
    
      delta = (E_StrBnd() - estbn);
      err = delta / estbn * 100;
      sprintf(logbuf, "Stretch-Bending        %8.3f      %8.3f     %8.3f     %5.2f", E_StrBnd(), estbn, delta, err);
      cout << logbuf << endl;
    
      delta = (E_OOP() - eoop);
      err = delta / eoop * 100;
      sprintf(logbuf, "Out-Of-Plane Bending   %8.3f      %8.3f     %8.3f     %5.2f", E_OOP(), eoop, delta, err);
      cout << logbuf << endl;
    
      delta = (E_Torsion() - etor);
      err = delta / etor * 100;
      sprintf(logbuf, "Torsional              %8.3f      %8.3f     %8.3f     %5.2f", E_Torsion(), etor, delta, err);
      cout << logbuf << endl;
    
      delta = (E_VDW() - evdw);
      err = delta / evdw * 100;
      sprintf(logbuf, "Van der Waals          %8.3f      %8.3f     %8.3f     %5.2f", E_VDW(), evdw, delta, err);
      cout << logbuf << endl;
      
      delta = (E_Electrostatic() - eeq);
      err = delta / ebond * 100;
      sprintf(logbuf, "Electrostatic          %8.3f      %8.3f     %8.3f     %5.2f", E_Electrostatic(), eeq, delta, err);
      cout << logbuf << endl;

      cout << endl;
      delta = (Energy() - etot);
      err = delta / etot * 100;
      sprintf(logbuf, "Total ENERGY           %8.3f      %8.3f     %8.3f     %5.2f", Energy(), etot, delta, err);
      cout << logbuf << endl;

    } // for (unsigned int c;; c++ )
   
    if (ifs)
      ifs.close();
    if (ifs2)
      ifs2.close();
  }
  
  int OBForceFieldMMFF94::GetBondType(OBAtom* a, OBAtom* b)
  {
    if (_mol.GetBond(a,b)->IsSingle()) {
      if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())) && !HasAromSet(atoi(a->GetType())) && !HasAromSet(atoi(b->GetType())))
        return 1;

      if (!IsInSameRing(a, b) && HasAromSet(atoi(a->GetType())) && HasAromSet(atoi(b->GetType())))
        return 1;
    } else {
/*
      if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())) && 
          a->IsAromatic() && b->IsAromatic() && 
	  !HasAromSet(atoi(a->GetType())) && !HasAromSet(atoi(b->GetType())) &&
	  IsInSameRing(a,b))
        return 1;
	*/
    }

    return 0;
  }
  
  int OBForceFieldMMFF94::GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int sumbondtypes;

    sumbondtypes = GetBondType(a,b) + GetBondType(b, c);

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3)) {
      switch (sumbondtypes) {
        case 0:
	  return 3; 
        case 1:
	  return 5; 
        case 2:
	  return 6; 
      }
    }
    
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4)) {
      switch (sumbondtypes) {
        case 0:
	  return 4; 
        case 1:
	  return 7; 
        case 2:
	  return 8; 
      }
    }

    return sumbondtypes;
  }
  
  int OBForceFieldMMFF94::GetStrBndType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int btab, btbc, atabc;
    bool inverse;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    atabc = GetAngleType(a, b, c);

    if (atoi(a->GetType()) <= atoi(c->GetType()))
      inverse = false;
    else
      inverse = true;

    switch (atabc) {
      case 0:
        return 0;

      case 1:
        if (btab)
	  if (!inverse)
	    return 1;
	  else
	    return 2;
	if (btbc)
	  if (!inverse)
	    return 2;
	  else
	    return 1;

      case 2:
        return 3;

      case 3:
        return 5;

      case 4:
        return 4;

      case 5:
        if (btab)
	  if (!inverse)
	    return 6;
	  else
	    return 7;
	if (btbc)
	  if (!inverse)
	    return 7;
	  else
	    return 6;
      
      case 6:
        return 8;
      
      case 7:
         if (btab)
	  if (!inverse)
	    return 9;
	  else
	    return 10;
	if (btbc)
	  if (!inverse)
	    return 10;
	  else
	    return 9;
      
      case 8:
        return 11;
    }
  }
  
  int OBForceFieldMMFF94::GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d)
  {
    int btab, btbc, btcd;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    btcd = GetBondType(c, d);
    
    if (btbc == 1)
      return 1;

    if (btbc && (btab || btcd))
      return 2;

    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && d->IsInRingSize(4))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d))
        return 4;

    if (a->IsInRingSize(5) && b->IsInRingSize(5) && c->IsInRingSize(5) && d->IsInRingSize(5))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d) && !a->IsAromatic())
        return 5; // needs checking for saturations

    return 0;
  }

  bool OBForceFieldMMFF94::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = _mol.GetSSSR();

    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    for (i = vr.begin();i != vr.end();i++) {
        for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
	    if (*j == a->GetIdx())
	        a_in = true;
	    if (*j == b->GetIdx())
	        b_in = true;
	}
        if (a_in && b_in)
            return true;
	a_in = false;
	b_in = false;
    }

    return false;
  }
  
  bool OBForceFieldMMFF94::HasLinSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _lin.begin(); i != _lin.end(); i++)
      if ((*i) == atomtype)
        return true;

    return false;
  }
 
  
  bool OBForceFieldMMFF94::HasAromSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _arom.begin(); i != _arom.end(); i++)
      if ((*i) == atomtype)
        return true;

    return false;
  }
 
  bool OBForceFieldMMFF94::HasSbmbSet(int atomtype)
  {
    vector<int>::iterator i;

    for (i = _sbmb.begin(); i != _sbmb.end(); i++)
      if ((*i) == atomtype)
        return true;

    return false;
  }

  OBFFParameter* OBForceFieldMMFF94::GetParameterMMFF94(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    if (!b)
      for (int idx=0; idx < parameter.size(); idx++)
        if ((a == parameter[idx].a) && (ffclass == parameter[idx].ipar5)) {
	  par = &parameter[idx];
	  return par;
	}

    if (!c)
      for (int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (ffclass == parameter[idx].ipar5)) || 
	    ((a == parameter[idx].b) && (b == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
	{
	  par = &parameter[idx];
	  return par;
	}

    if (!d)
      for (int idx=0; idx < parameter.size(); idx++)
        if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (ffclass == parameter[idx].ipar5)) || 
	    ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
	{
	  par = &parameter[idx];
	  return par;
        }

    for (int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (d == parameter[idx].d) && (ffclass == parameter[idx].ipar5)) || 
          ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && (d == parameter[idx].a) && (ffclass == parameter[idx].ipar5))) 
      {
	par = &parameter[idx];
	return par;
      }

    return NULL;
  }
 
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field