/**********************************************************************
Copyright (C) 2005,2006,2007 Chris Morley

Based on the IUPAC InChI reference software, which is distributed
under the GNU LGPL:
Copyright (C) 2005 The International Union of Pure and Applied Chemistry
IUPAC International Chemical Identifier (InChI) (contact:secretariat@iupac.org)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include "inchi_api.h"
#include <sstream>
#include <set>
#include <vector>
#include <iterator>
#include <openbabel/inchiformat.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>

using namespace std;
namespace OpenBabel
{
extern string GetInChI(istream& is);

//Make an instance of the format class
InChIFormat theInChIFormat;

/////////////////////////////////////////////////////////////////
bool InChIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL) return false;
  istream &ifs = *pConv->GetInStream();

  //Extract InChI from input stream even if it is split
  string inchi;
  do
  {
    inchi = GetInChI(ifs);
    if(inchi.empty())
      return false; //eof
  }while(inchi.size()<9); //ignore empty "InChI=" or "InChI=1/"

  // Save the original inchi to be used on output,
  // avoiding conversion to and from OBMol
  SaveInchi(pmol, inchi);

  //Set up input struct
  inchi_InputINCHI inp;

  char* opts= GetInChIOptions(pConv, true);
  inp.szOptions = opts;

  char* nonconstinchi = new char[inchi.size()+1];
  inp.szInChI = strcpy(nonconstinchi, inchi.c_str());

  inchi_OutputStruct out;
  memset(&out, 0, sizeof(out));

  //Call the conversion routine in InChI code
  int ret = GetStructFromINCHI( &inp, &out );

  if (ret!=inchi_Ret_OKAY)
  {
    string mes = out.szMessage ? out.szMessage : "Conversion failed";
    obErrorLog.ThrowError("InChI code", mes + '\n' + inchi, obWarning);
  }
  delete[] nonconstinchi;
  delete[] opts;
  if(ret)
    return false;

  //Read name if requested e.g InChI=1/CH4/h1H4 methane
  //OR InChI=1/CH4/h1H4 "First alkane"  Quote can be any punct char and
  //uses upto the end of the line if second quote is not found
  if(pConv->IsOption("n",OBConversion::INOPTIONS))
  {
    string name;
    if(getline(ifs, name))
      pmol->SetTitle(Trim(name));
  }

  //Option to add InChI text to title
  if(pConv->IsOption("a",OBConversion::INOPTIONS))
  {
    string title(pmol->GetTitle());
    title += ' ' + inchi;
    pmol->SetTitle(title);
  }

  //Translate the returned structure into OBMol
  pmol->SetDimension(0);
  pmol->BeginModify();
  int i;
  //Make all atoms first because we need pointers to later ones
  for(i=0;i<out.num_atoms;++i)
  {
    pmol->NewAtom();
  }
  for(i=0;i<out.num_atoms;++i)
  {
    OBAtom* patom = pmol->GetAtom(i+1); //index starts at 1
    inchi_Atom* piat = &out.atom[i];
    int iso=0;
    patom->SetAtomicNum(etab.GetAtomicNum(piat->elname,iso));
    patom->SetIsotope(iso);
    if(piat->isotopic_mass)
      patom->SetIsotope(piat->isotopic_mass - ISOTOPIC_SHIFT_FLAG +
          (int)(isotab.GetExactMass(patom->GetAtomicNum())+0.5));

    patom->SetSpinMultiplicity(piat->radical);
    patom->SetFormalCharge(piat->charge);
//    patom->SetVector(piat->x,piat->y,piat->z);

    int j;
    for(j=0;j<piat->num_bonds;++j)
    {
      pmol->AddBond(i+1, piat->neighbor[j]+1, piat->bond_type[j]);
    }

    //OB takes care of implicit hydrogens, so num_iso_H[0] and num_iso_H[1] are ignored,
    //except for H2, which has one explicit H and one implict H. OB doesn't add implicit H to H.
    //Implicit D and T also need to be added explicitly.
    for(int m=0;m<3;++m)
    {
      if(piat->num_iso_H[m] && (m>1 || *piat->elname=='H'))
      {
        for(int k=0;k<piat->num_iso_H[m];++k)
        {
          OBAtom* DorT = pmol->NewAtom();
          DorT->SetAtomicNum(1);
          DorT->SetIsotope(m);
          pmol->AddBond(i+1, pmol->NumAtoms(), 1);
        }
      }
    }
  }

  //***@todo implicit H isotopes
  //Stereochemistry
  for(i=0;i<out.num_stereo0D;++i)
  {
    inchi_Stereo0D& stereo = out.stereo0D[i];

    switch(stereo.type)
    {
    case INCHI_StereoType_DoubleBond:
    {
      OBCisTransStereo::Config *ct = new OBCisTransStereo::Config;

      ct->begin = stereo.neighbor[1];
      ct->end = stereo.neighbor[2];
      unsigned long start = OBStereo::ImplicitRef;
      unsigned long end = OBStereo::ImplicitRef;
      FOR_NBORS_OF_ATOM(a, pmol->GetAtom(ct->begin + 1)) {
        if ( !(a->GetId() == ct->end || a->GetId() == stereo.neighbor[0] ) ) {
          start = a->GetId();
          break;
        }
      }
      FOR_NBORS_OF_ATOM(b, pmol->GetAtom(ct->end + 1)) {
        if ( !(b->GetId() == ct->begin || b->GetId() == stereo.neighbor[3] ) ) {
          end = b->GetId();
          break;
        }
      }
      ct->refs = OBStereo::MakeRefs(stereo.neighbor[0], start, stereo.neighbor[3], end);

      if(stereo.parity==INCHI_PARITY_EVEN)
        ct->shape = OBStereo::ShapeU;
      else if(stereo.parity==INCHI_PARITY_ODD)
        ct->shape = OBStereo::ShapeZ;
      else
        ct->specified = false;

      OBCisTransStereo *obct = new OBCisTransStereo(pmol);
      obct->SetConfig(*ct);
      pmol->SetData(obct);

      //InChI seems to prefer to define double bond stereo using H atoms.

      break;
    }
    case INCHI_StereoType_Tetrahedral:
    {
      OBTetrahedralStereo::Config *ts = new OBTetrahedralStereo::Config;
      ts->center = stereo.central_atom;
      ts->from = stereo.neighbor[0];
      if (ts->from == ts->center) // Handle the case where there are only three neighbours
        ts->from = OBStereo::ImplicitRef;
      ts->refs = OBStereo::MakeRefs(stereo.neighbor[1], stereo.neighbor[2],
                                    stereo.neighbor[3]);

      if(stereo.parity==INCHI_PARITY_EVEN)
        ts->winding = OBStereo::Clockwise;
      else if(stereo.parity==INCHI_PARITY_ODD)
        ts->winding = OBStereo::AntiClockwise;
      else
        ts->specified = false;

      OBTetrahedralStereo *obts = new OBTetrahedralStereo(pmol);
      obts->SetConfig(*ts);
      pmol->SetData(obts);

      break;
    }

    case INCHI_StereoType_Allene:
    default:
      obErrorLog.ThrowError("InChI code", "Unsupported stereo type has been ignored.", obWarning);
    }
  }
  // Tidy up the stereo chemistry by removing any objects that are not
  // consistent with OB's symmetry analysis
  StereoFrom0D(pmol);

  pmol->EndModify();
  FreeStructFromINCHI( &out );
  return true;
}

int InChIFormat::SkipObjects(int n, OBConversion* pConv)
{
  istream& ifs = *pConv->GetInStream();
  string inchi;
  while(ifs.good() && n)
  {
    inchi = GetInChI(ifs);
    if(inchi.size()>=8)//ignore empty "InChI=" or "InChI=1/"
      --n;
  }
  return ifs.good() ? 1 : -1;
}

/////////////////////////////////////////////////////////////////
bool InChIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  //Although the OBMol may be altered, it is restored before exit.
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL) return false;
    OBMol& mol = *pmol;

  string ostring; //the inchi string
  inchi_Output inout;

  stringstream molID;
  if(strlen(mol.GetTitle())==0)
    molID << '#' << pConv->GetOutputIndex() << ' ';
  else
    molID << mol.GetTitle() << ' ';
  if(pConv->GetOutputIndex()==1)
    firstID=molID.str();

  //Use any existing InChI, probably from an InChIformat input,
  //in preference to determining it from the structure.
  if(!pConv->IsOption("r") && pmol->HasData("inchi"))
  {
    //All origins for the data are currently acceptable.
    //Possibly this may need to be restricted to data with a local origin.
    ostring = pmol->GetData("inchi")->GetValue();
  }
  else
  {
    //Determine InChI from the chemical structure
    if(pmol->NumAtoms()==0) return true; // PR#2864334

    inchi_Input inp;
    memset(&inp,0,sizeof(inchi_Input));

    // Prepare stereo information for 2D, 3D
    map<OBBond*, OBStereo::BondDirection> updown;
    map<OBBond*, OBStereo::Ref> from;
    map<OBBond*, OBStereo::Ref>::const_iterator from_cit;
    if (mol.GetDimension() == 3 || (mol.GetDimension()==2 && pConv->IsOption("s", pConv->OUTOPTIONS)!=NULL))
      TetStereoToWedgeHash(mol, updown, from);
    set<OBBond*> unspec_ctstereo = GetUnspecifiedCisTrans(mol);

    OBAtom* patom;
    vector<inchi_Atom> inchiAtoms(mol.NumAtoms());
    vector<OBNodeBase*>::iterator itr;
    for (patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
    {
      //OB atom index starts at 1; inchi atom index starts at 0
      inchi_Atom& iat = inchiAtoms[patom->GetIdx()-1];
      memset(&iat,0,sizeof(inchi_Atom));

      iat.x = patom->GetX();
      iat.y = patom->GetY();
      iat.z = patom->GetZ();

      int nbonds = 0;
      vector<OBBond*>::iterator itr;
      OBBond *pbond;
      for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr))
      {
        from_cit = from.find(pbond);
        // Do each bond only once. If the bond is a stereobond
        // ensure that the BeginAtom is the 'from' atom.
        if( (from_cit==from.end() && patom->GetIdx() != pbond->GetBeginAtomIdx()) ||
            (from_cit!=from.end() && from_cit->second != patom->GetId()) )
          continue;

        iat.neighbor[nbonds]      = pbond->GetNbrAtomIdx(patom)-1;
        int bo = pbond->GetBO();
        if(bo==5)
          bo=4;
        iat.bond_type[nbonds]     = bo;

        OBStereo::BondDirection stereo = OBStereo::NotStereo;
        if (mol.GetDimension()==2 && pConv->IsOption("s", pConv->OUTOPTIONS)==NULL) {
          if (pbond->IsWedge())
            stereo = OBStereo::UpBond;
          else if (pbond->IsHash())
            stereo = OBStereo::DownBond;
          else if (pbond->IsWedgeOrHash())
            stereo = OBStereo::UnknownDir;
        } 
        else if (from_cit!=from.end()) { // It's a stereo bond
          stereo = updown[pbond];
        }

        if (mol.GetDimension() != 0) {
          inchi_BondStereo2D bondstereo2D = INCHI_BOND_STEREO_NONE;
          if (stereo != OBStereo::NotStereo) {
            switch (stereo) {
              case OBStereo::UpBond:
                bondstereo2D = INCHI_BOND_STEREO_SINGLE_1UP;
                break;
              case OBStereo::DownBond:
                bondstereo2D = INCHI_BOND_STEREO_SINGLE_1DOWN;
                break;
              case OBStereo::UnknownDir:
                bondstereo2D = INCHI_BOND_STEREO_SINGLE_1EITHER;
                break;
              default:
                ; // INCHI_BOND_STEREO_NONE
            }
          }
          // Is it a double bond with unspecified stereochemistry?
          if (unspec_ctstereo.find(pbond)!=unspec_ctstereo.end())
            bondstereo2D = INCHI_BOND_STEREO_DOUBLE_EITHER;
          iat.bond_stereo[nbonds] = bondstereo2D;
        }
        nbonds++;
      }

      strcpy(iat.elname,etab.GetSymbol(patom->GetAtomicNum()));
      iat.num_bonds = nbonds;
      //Let inchi add implicit Hs unless the atom is known not to have any
      iat.num_iso_H[0] = patom->HasNoHForced() || patom->IsMetal() ? 0 : -1;
      if(patom->GetIsotope())
      {
        iat.isotopic_mass = ISOTOPIC_SHIFT_FLAG +
          patom->GetIsotope() - (int)(etab.GetMass(patom->GetAtomicNum())+0.5);
      }
      else
        iat.isotopic_mass = 0 ;
      iat.radical = patom->GetSpinMultiplicity();
      //InChI doesn't recognize spin miltiplicity of 4 or 5 (as used in OB for CH and C atom)
      if(iat.radical>=4)
      {
        iat.radical=0;
        iat.num_iso_H[0] = 0; //no implicit hydrogens
      }
      iat.charge  = patom->GetFormalCharge();
    }

    inp.atom = &inchiAtoms[0];

    vector<inchi_Stereo0D> stereoVec;

    if(mol.GetDimension()==0)
    {
      std::vector<OBGenericData*>::iterator data;
      std::vector<OBGenericData*> stereoData = mol.GetAllData(OBGenericDataType::StereoData);
      for (data = stereoData.begin(); data != stereoData.end(); ++data) {
        if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::Tetrahedral) {
          OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
          OBTetrahedralStereo::Config config = ts->GetConfig();

          if(config.specified) {
            inchi_Stereo0D stereo;
            stereo.type = INCHI_StereoType_Tetrahedral;
            stereo.central_atom = static_cast<AT_NUM> (config.center);

            bool has_implicit = false; // Does chirality involve implicit lone pair?
            if (config.from == OBStereo::ImplicitRef || config.refs[0] == OBStereo::ImplicitRef) {
              has_implicit = true;
              config = ts->GetConfig(OBStereo::ImplicitRef); // Make the 'from' atom the lone pair
            }

            if (!has_implicit)
              stereo.neighbor[0] = static_cast<AT_NUM> (config.from);
            else
              stereo.neighbor[0] = stereo.central_atom;
            for(int i=0; i<3; ++i)
              stereo.neighbor[i + 1] = static_cast<AT_NUM> (config.refs[i]);

            if (config.winding == OBStereo::Clockwise)
              stereo.parity = INCHI_PARITY_EVEN;
            else
              stereo.parity = INCHI_PARITY_ODD;

            stereoVec.push_back(stereo);
          }
        }
      }

      //Double bond stereo (still inside 0D section)
      //Currently does not handle cumulenes
      for (data = stereoData.begin(); data != stereoData.end(); ++data) {
        if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::CisTrans) {
          OBCisTransStereo *ts = dynamic_cast<OBCisTransStereo*>(*data);
          OBCisTransStereo::Config config = ts->GetConfig();

          if(config.specified) {
            inchi_Stereo0D stereo;
            stereo.central_atom = NO_ATOM;
            stereo.type = INCHI_StereoType_DoubleBond;
            OBStereo::Refs refs = config.refs;
            unsigned long start = refs[0];
            if (refs[0]==OBStereo::ImplicitRef)
              start = refs[1];
            unsigned long end = refs[3];
            if (refs[3]==OBStereo::ImplicitRef)
              end = refs[2];

            stereo.neighbor[0] = static_cast<AT_NUM> (start);
            stereo.neighbor[1] = static_cast<AT_NUM> (config.begin);
            stereo.neighbor[2] = static_cast<AT_NUM> (config.end);
            stereo.neighbor[3] = static_cast<AT_NUM> (end);

            if (ts->IsTrans(start, end))
              stereo.parity = INCHI_PARITY_EVEN;
            else
              stereo.parity = INCHI_PARITY_ODD;

            stereoVec.push_back(stereo);
          }
        }
      }
    }

    char* opts = GetInChIOptions(pConv, false);
    inp.szOptions = opts;

    inp.num_atoms = mol.NumAtoms();
    inp.num_stereo0D = stereoVec.size();
    if(inp.num_stereo0D>0)
      inp.stereo0D = &stereoVec[0];

    //inchi_Output inout; now declared in block above
    memset(&inout,0,sizeof(inchi_Output));

    int ret = GetINCHI(&inp, &inout);

    delete[] opts;
    if(ret!=inchi_Ret_OKAY)
    {
      if(inout.szMessage)
      {
        string mes(inout.szMessage);
        if(pConv->IsOption("w"))
        {
          string::size_type pos;
          string targ[4];
          targ[0] = "Omitted undefined stereo";
          targ[1] = "Charges were rearranged";
          targ[2] = "Proton(s) added/removed";
          targ[3] = "Metal was disconnected";
          for(int i=0;i<4;++i)
          {
            pos = mes.find(targ[i]);
            if(pos!=string::npos)
            {
              mes.erase(pos,targ[i].size());
              if(mes[pos]==';')
                mes[pos]=' ';
            }
          }
        }
        Trim(mes);
        if(!mes.empty())
          obErrorLog.ThrowError("InChI code", molID.str() + ':' + mes, obWarning);
      }

      if(ret!=inchi_Ret_WARNING)
      {
        obErrorLog.ThrowError("InChI code", "InChI generation failed", obError);
        FreeStdINCHI(&inout);
        return false;
      }
    }
    ostring = inout.szInChI;
  }

  //Truncate the InChI if requested
  const char* truncspec = pConv->IsOption("T");
  if(truncspec)
  {
    string trunc(truncspec);
    EditInchi(ostring, trunc);
  }

  if(pConv->IsOption("K")) //Generate InChIKey and add after InChI on same line
  {
    char szINCHIKey[28];
    GetStdINCHIKeyFromStdINCHI(ostring.c_str(), szINCHIKey);
    ostring = szINCHIKey;
  }

  if(pConv->IsOption("t"))
  {
    ostring += ' ';
    ostring +=  pmol->GetTitle();
  }

  ostream &ofs = *pConv->GetOutStream();

  if(pConv->IsOption("U"))
  {
    if(pConv->GetOutputIndex()==1)
      allInchi.clear();
    //Just add to set and don't output, except at the end
    allInchi.insert(ostring);

    if(pConv->IsLast())
    {
      nSet::iterator itr;
      for(itr=allInchi.begin();itr!=allInchi.end();++itr)
        ofs << *itr << endl;
    }
    return true;
  }
  else if(pConv->IsOption("u"))
  {
    if(pConv->GetOutputIndex()==1)
      allInchi.clear();
    if(!allInchi.insert(ostring).second)
      return true; //no output if already in set
  }

  ofs << ostring << endl;

  if (pConv->IsOption("a"))
    ofs << inout.szAuxInfo << endl;

  if(pConv->IsOption("l"))
    //Display InChI log message. With multiple molecules, it appears only once
    obErrorLog.ThrowError("InChI log", inout.szLog, obError, onceOnly);

  if(pConv->IsOption("e"))
  {
    if(pConv->GetOutputIndex()==1)
      firstInchi = inout.szInChI;
    else
    {
      ofs << "Molecules " << firstID << "and " << molID.str();
      ofs << InChIErrorMessage(CompareInchi(firstInchi.c_str(), inout.szInChI)) << endl;
    }
  }

  FreeStdINCHI(&inout);
  return true;
}

//////////////////////////////////////////////////////////
char InChIFormat::CompareInchi(const string& Inchi1, const string& Inchi2)
{
  //Returns 0 if identical or an char identifying the layer where they first differed
  string s1(Inchi1), s2(Inchi2);

  if(s1.size()<s2.size())
    s1.swap(s2);
  string::size_type pos;
  for(pos=0;pos<s1.size();++pos)
  {
    if(pos==s2.size() || s1[pos]!=s2[pos])
      return s1[s1.rfind('/',pos)+1];
  }
  return 0;

/*  //Remove anything after the end of the Inchi
  string::size_type pos;
  pos = s1.find_first_of(" \t\n");
  if(pos!=string::npos)
    s1.erase(pos);
  pos = s2.find_first_of(" \t\n");
  if(pos!=string::npos)
    s2.erase(pos);

  vector<string> layers1, layers2;
  tokenize(layers1,s1,"/\n");
  tokenize(layers2,s2,"/\n");
  unsigned int i;
  if(layers1.size()<layers2.size())
    layers1.swap(layers2); //layers1 is the longest

  for(i=1;i<layers2.size();++i)
  {
    if(layers1[i]!=layers2[i])
    {
      char ch = '+';
      if(i>1) //not formula layer
        ch=layers1[i][0];
      return ch;
    }
  }
  if(layers1.size()==layers2.size())
    return 0;
  else
    return layers1[i][0];
*/
}

string InChIFormat::InChIErrorMessage(const char ch)
{
  string s;
  switch (ch)
  {
  case 0:
    s = " are identical";
    break;
  case '+':
    s = " have different formulae";
    break;
  case 'c':
    s = " have different connection tables";
    break;
  case 'h':
    s = " have different bond orders, or radical character";
    break;
  case 'q':
    s = " have different charges";
    break;
  case 'p':
    s = " have different numbers of attached protons";
    break;
  case 'b':
    s = " have different double bond stereochemistry";
    break;
  case 'm':
  case 't':
    s = " have different sp3 stereochemistry";
    break;
  case 'i':
    s = " have different isotopic composition";
    break;
  default:
    s = " are different";
  }
  return s;
}

OBAtom* InChIFormat::GetCommonAtom(OBBond* pb1, OBBond* pb2)
{
  OBAtom* pa1 = pb1->GetBeginAtom();
  if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
    return pa1;
  pa1 = pb1->GetEndAtom();
  if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
    return pa1;
  return NULL; //not adjacent bonds
}


//Returns pointer to InChI options string, which needs to be deleted with delete[]
//If there are no options returns an empty string
char* InChIFormat::GetInChIOptions(OBConversion* pConv, bool Reading)
{
  vector<string> optsvec; //the InChi options
  /* In Standard InChI these are not used
  if(!Reading && !pConv->IsOption("n"))
    //without -xn option, the default is to write using these 'recommended' options
    tokenize(optsvec, "FixedH RecMet SPXYZ SAsXYZ Newps Fb Fnud");
  */
  char* opts;
  OBConversion::Option_type opttyp = Reading ? OBConversion::INOPTIONS : OBConversion::OUTOPTIONS;
  const char* copts = pConv->IsOption("X", opttyp);
  if(copts)
  {
    string tmp(copts); // GCC doesn't like passing string temporaries to functions
    vector<string> useropts;
    tokenize(useropts, tmp);
    copy(useropts.begin(), useropts.end(), back_inserter(optsvec));
  }

  //Add a couple InChI options built in to OB
  if(opttyp==OBConversion::OUTOPTIONS)
  {
    if(pConv->IsOption("F", opttyp))
    {
      string tmp2("FixedH");
      optsvec.push_back(tmp2);
    }
    if(pConv->IsOption("M", opttyp))
    {
      string tmp2("RecMet");
      optsvec.push_back(tmp2);
    }
  }


#ifdef WIN32
    string ch(" /");
#else
    string ch(" -");
#endif

  string sopts;
  for(int i=0;i<optsvec.size();++i)
    sopts += ch + optsvec[i];
  opts = new char[strlen(sopts.c_str())+1]; //has to be char, not const char
  return strcpy(opts, sopts.c_str());
  opts = new char[1];
  *opts = '\0';
  return opts;
}

bool InChIFormat::EditInchi(std::string& inchi, std::string& spec)
{
  std::vector<std::string> vec;
  std::vector<std::string>::iterator itr;
  tokenize(vec, spec, " \t/");
  for(itr=vec.begin();itr!=vec.end();++itr)
  {
    if(*itr=="formula")
    {
      std::string::size_type pos = inchi.find('/', inchi.find('/')+1); //2nd /
      if(pos!=string::npos)
        inchi.erase(pos);
    }
    else if(*itr=="connect")
      RemoveLayer(inchi,"/h",true);
    else if(*itr=="nochg")
    {
      RemoveLayer(inchi,"/p");
      RemoveLayer(inchi,"/q");
    }
    else if(*itr=="nosp3")
    {
      RemoveLayer(inchi,"/t");
      RemoveLayer(inchi,"/m");
      RemoveLayer(inchi,"/s");
    }
    else if(*itr=="noEZ")
      RemoveLayer(inchi,"/b");
    else if(*itr=="noiso")
      RemoveLayer(inchi,"/i");
    else if(*itr=="nostereo")
    {
      RemoveLayer(inchi,"/t");
      RemoveLayer(inchi,"/m");
      RemoveLayer(inchi,"/s");
      RemoveLayer(inchi,"/b");
    }
    else if(!(*itr).empty())
    {
      obErrorLog.ThrowError(__FUNCTION__,
      spec + " not recognized as a truncation specification",obError, onceOnly);
      return false;
    }
  }
  return true;
}

void InChIFormat::RemoveLayer (std::string& inchi, const std::string& str, bool all)
{
  std::string::size_type pos = inchi.find(str);
  if(pos!=string::npos)
    inchi.erase(pos, (all ? string::npos : inchi.find('/', pos+1) - pos));
}

void InChIFormat::SaveInchi(OBMol* pmol, const std::string& s)
{
  OBPairData* dp = new OBPairData;
  dp->SetAttribute("inchi");
  dp->SetValue(s);
  dp->SetOrigin(local);
  pmol->SetData(dp);
}

//************************************************************************
InChICompareFormat theInChICompareFormat;

bool InChICompareFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  pConv->AddOption("e",OBConversion::OUTOPTIONS);
  pConv->AddOption("t",OBConversion::OUTOPTIONS);
  return theInChIFormat.WriteMolecule(pOb,pConv);
}

}//namespace OpenBabel

