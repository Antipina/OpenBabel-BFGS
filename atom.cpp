/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
	  
#include "mol.h"
#include "typer.h"
#include "molchrg.h"
#include "phmodel.h"

namespace OpenEye {

extern OEAromaticTyper  aromtyper;
extern OEAtomTyper      atomtyper;
extern OEPhModel        phmodel;

//
// OEAtom member functions
//

OEAtom::OEAtom()
{
  _c = (float**)NULL;
  _parent = (OEMol*)NULL;
  _cidx = 0;
  _flags=0;
  _idx = 0;
  _hyb = 0;
  _ele = (char)0;
  _impval = 0;
  _fcharge = 0;
  //_stereo = 0;
  _type[0] = '\0';
  _pcharge = 0.0;
  _vbond.clear();
  _vbond.reserve(4);
  _residue = (OEResidue*)NULL;
}

OEAtom::~OEAtom()
{
    if (_residue != NULL)
        _residue->RemoveAtom(this);
}

void OEAtom::Clear()
{
  _c = (float**)NULL;
  _cidx = 0;
  _flags=0;
  _idx = 0;
  _hyb = 0;
  _ele = (char)0;
  _impval = 0;
  _fcharge = 0;
  //_stereo = 0;
  _type[0] = '\0';
  _pcharge = 0.0;
  _vbond.clear();
  _vbond.reserve(4);
  _residue = (OEResidue*)NULL;
}

OEAtom &OEAtom::operator=(OEAtom &src)
     //copy atom information 
     //bond info is not copied here as ptrs may be invalid
{
  _idx = src.GetIdx();
  _hyb = src.GetHyb();
  _ele = src.GetAtomicNum();
  _fcharge = src.GetFormalCharge();
  //_stereo = src.GetStereo();
  strcpy(_type,src.GetType());
  _pcharge = src.GetPartialCharge();
  _v = src.GetVector();
  _flags = src.GetFlag();
  _residue = (OEResidue*)NULL;
  return(*this);
}

bool OEAtom::IsConnected(OEAtom *a1)
{
  vector<OEBond*>::iterator i;
  OEBond *bond;

  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBeginAtom() == a1 || bond->GetEndAtom() == a1)
      return(true);

  return(false);
}

bool OEAtom::IsOneThree(OEAtom *a1)
{
  OEAtom *atom1,*atom2;
  OEBond *bond1,*bond2;
  vector<OEBond*>::iterator i,j;
  atom1 = this;
  atom2 = a1;

  for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
    for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
      if (bond1->GetNbrAtom(atom1) == bond2->GetNbrAtom(atom2))
	return(true);

  return(false);
}

bool OEAtom::IsOneFour(OEAtom *a1)
{
  OEAtom *atom1,*atom2;
  OEBond *bond1,*bond2;
  vector<OEBond*>::iterator i,j;
  atom1 = this;
  atom2 = a1;

  for (bond1 = atom1->BeginBond(i);bond1;bond1 = atom1->NextBond(i))
    for (bond2 = atom2->BeginBond(j);bond2;bond2 = atom2->NextBond(j))
      if ((bond1->GetNbrAtom(atom1))->IsConnected(bond2->GetNbrAtom(atom2)))
	return(true);

  return(false);
}

bool OEAtom::IsAxial()
{
  float tor;
  OEAtom *a,*b,*c;
  vector<OEBond*>::iterator i,j,k;
  
  for (a = BeginNbrAtom(i);a;a = NextNbrAtom(i))
    if (a->GetHyb() == 3 && a->IsInRing() && !(*i)->IsInRing())
      for (b = a->BeginNbrAtom(j);b;b = a->NextNbrAtom(j))
	if (b != this && b->IsInRing() && b->GetHyb() == 3)
	  for (c = b->BeginNbrAtom(k);c;c = b->NextNbrAtom(k))
	    if (c != a && c->IsInRing())
	      {
			tor = fabs(((OEMol*)GetParent())->GetTorsion(this,a,b,c));
			return(tor > 55.0 && tor < 75.0);
	      }

  return(false);
}


bool OEAtom::HasAlphaBetaUnsat(bool includePandS)
{
  OEAtom *a1,*a2;
  vector<OEBond*>::iterator i,j;

  for (a1 = BeginNbrAtom(i);a1;a1 = NextNbrAtom(i))
    if (includePandS || (!a1->IsPhosphorus() && !a1->IsSulfur()))
    for (a2 = a1->BeginNbrAtom(j);a2;a2 = a1->NextNbrAtom(j))
      if (a2 != this && ((*j)->GetBO() == 2 || (*j)->GetBO() == 3 || (*j)->GetBO() == 5))
	return(true);

  return(false);
}

bool OEAtom::HasBondOfOrder(int order)
{
  OEBond *bond;
  vector<OEBond*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBO() == order)
      return(true);

  return(false);
}

int OEAtom::CountBondsOfOrder(int order)
{
	int count = 0;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if (bond->GetBO() == order)
      count++;

  return(count);
}

bool OEAtom::IsPolarHydrogen()
{
  if (!IsHydrogen()) return(false);
  
  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      atom = bond->GetNbrAtom(this);
      if (atom->GetAtomicNum() == 7) return(true);
      if (atom->GetAtomicNum() == 8) return(true);
      if (atom->GetAtomicNum() == 15) return(true);
      if (atom->GetAtomicNum() == 16) return(true);
    }

  return(false);
}

bool OEAtom::IsNonPolarHydrogen()
{
  if (!IsHydrogen()) return(false);
  
  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      atom = bond->GetNbrAtom(this);
      if (atom->GetAtomicNum() == 6) return(true);
    }

  return(false);
}

Vector &OEAtom::GetVector()
{
  if (!_c) return(_v);

  _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
  return(_v);
}

void OEAtom::SetVector()
{
  oeAssert(_c);
  if (_c) _v.Set((*_c)[_cidx],(*_c)[_cidx+1],(*_c)[_cidx+2]);
}

void OEAtom::SetVector(Vector &v)
{
  if (!_c) _v = v;
  else
    {
      (*_c)[_cidx  ] = v.x();
      (*_c)[_cidx+1] = v.y();
      (*_c)[_cidx+2] = v.z();
    }
}

void OEAtom::SetVector(const float x,const float y,const float z)
{
  if (!_c) _v.Set(x,y,z);
  else
    {
      (*_c)[_cidx  ] = x;
      (*_c)[_cidx+1] = y;
      (*_c)[_cidx+2] = z;
    }
}

OEAtom *OEAtom::GetNextAtom()
{
  OEMol *mol = (OEMol*)GetParent();
  return(((unsigned)GetIdx() == mol->NumAtoms())? NULL : mol->GetAtom(GetIdx()+1));
}

OEResidue *OEAtom::GetResidue()
{
    if (_residue != NULL)
        return _residue;
    else if (!((OEMol*)GetParent())->HasChainsPerceived())
    {
        chainsparser.PerceiveChains(*((OEMol*)GetParent()));
        return _residue;
    }
    else
        return NULL;
}

char *OEAtom::GetType()
{
  OEMol *mol = (OEMol*)GetParent();
  if (mol && !mol->HasAtomTypesPerceived())
    atomtyper.AssignTypes(*((OEMol*)GetParent()));

  return(_type);
}

unsigned int OEAtom::GetImplicitValence() const
{
  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();
  if (mol && !mol->HasImplicitValencePerceived())
    atomtyper.AssignImplicitValence(*((OEMol*)((OEAtom*)this)->GetParent()));

  return((unsigned int)_impval);
}

unsigned int OEAtom::GetHyb() const
{
  //hybridization is assigned when atoms are typed
  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();
  if (mol && !mol->HasHybridizationPerceived())
    atomtyper.AssignHyb(*mol);

  return(_hyb);
}


unsigned int OEAtom::GetHvyValence() const
     //returns the number of non-hydrogens connected to an atom
{
  unsigned int count=0;

  OEAtom *atom;
  vector<OEBond*>::iterator i;
  for (atom = ((OEAtom*)this)->BeginNbrAtom(i);atom;atom = ((OEAtom*)this)->NextNbrAtom(i))
	  if (!atom->IsHydrogen())
		  count++;

  return(count);
}

unsigned int OEAtom::GetHeteroValence() const
     //returns the number of heteroatoms connected to an atom
{
  unsigned int count=0;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  for (bond = ((OEAtom*)this)->BeginBond(i);bond;bond = ((OEAtom*)this)->NextBond(i))
    if (bond->GetNbrAtom((OEAtom*)this)->IsHeteroatom())
      count++;

  return((unsigned int)count);
}

float OEAtom::GetPartialCharge()
{
  if (!GetParent()) return(_pcharge);
  if (!((OEMol*)GetParent())->AutomaticPartialCharge()) return(_pcharge);

  if (!((OEMol*)GetParent())->HasPartialChargesPerceived())
    {
      //seed partial charges are set in the atom typing procedure
      phmodel.AssignSeedPartialCharge(*((OEMol*)GetParent()));
      OEGastChrg gc;
      gc.AssignPartialCharges(*((OEMol*)GetParent()));
    }

  return(_pcharge);
}

bool OEAtom::IsAmideNitrogen()
     //returns true if nitrogen is part of an amide
{
  if (!IsNitrogen()) return(false);

  OEAtom *nbratom,*atom;
  OEBond *abbond,*bond;

  vector<OEBond*>::iterator i,j;
  atom = this;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    {
      nbratom = bond->GetNbrAtom(atom);
      for (abbond = nbratom->BeginBond(j);abbond;abbond = nbratom->NextBond(j))
	if (abbond->GetBO() == 2 && 
	   (((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 8) ||
			((abbond->GetNbrAtom(nbratom))->GetAtomicNum() == 16)))
	  return(true);
    }  

  return(false);
}

bool OEAtom::IsAromaticNOxide()
{
	if (!IsNitrogen() || !IsAromatic()) return(false);

	OEAtom *atom;
	vector<OEBond*>::iterator i;

	for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
		if (atom->IsOxygen() && !(*i)->IsInRing() && (*i)->GetBO() == 2)		
			return(true);

	return(false);
}

bool OEAtom::IsCarboxylOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsCarbon())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() != 2) return(false);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OEAtom::IsPhosphateOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsPhosphorus())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() > 2) return(true);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(false);
}

bool OEAtom::IsSulfateOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsSulfur())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() < 3) return(false);

  //atom is connected to a carbon that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OEAtom::IsNitroOxygen()
{
  if (!IsOxygen()) return(false);
  if (GetHvyValence() != 1) return(false);

  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;

  atom = NULL;
  for (bond = BeginBond(i);bond;bond = NextBond(i))
    if ((bond->GetNbrAtom(this))->IsNitrogen())
      {
	atom = bond->GetNbrAtom(this);
	break;
      }
  if (!atom) return(false);
  if (atom->CountFreeOxygens() != 2) return(false);

  //atom is connected to a nitrogen that has a total 
  //of 2 attached free oxygens
  return(true);
}

bool OEAtom::IsHeteroatom()
{
	switch(GetAtomicNum())
	{
	case 7:
	case 8:
	case 9:
	case 15:
	case 16:
	case 17:
	case 35:
	case 53:
		return(true);
	}
	return(false);
}

bool OEAtom::IsAromatic() const
{
  if (((OEAtom*)this)->HasFlag(OE_AROMATIC_ATOM)) return(true);

  OEMol	*mol = (OEMol*)((OEAtom*)this)->GetParent();

  if (!mol->HasAromaticPerceived())
    {
      aromtyper.AssignAromaticFlags(*mol);
      if (((OEAtom*)this)->HasFlag(OE_AROMATIC_ATOM)) return(true);
    }

  return(false);
}

bool OEAtom::IsInRing() const
{
  if (((OEAtom*)this)->HasFlag(OE_RING_ATOM)) return(true);

  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();
  if (!mol->HasRingAtomsAndBondsPerceived())
    {
      mol->FindRingAtomsAndBonds();
      if (((OEAtom*)this)->HasFlag(OE_RING_ATOM)) return(true);
    }

  return(false);
}

bool OEAtom::IsChiral()
{
  if (HasFlag(OE_CHIRAL_ATOM)) return(true);
  
  if (!((OEMol*)GetParent())->HasChiralityPerceived())
    {
      ((OEMol*)GetParent())->FindChiralCenters();
      if (HasFlag(OE_CHIRAL_ATOM)) return(true);
    }

  return(false);
}

bool OEAtom::IsInRingSize(int size) const
{
  vector<OERing*> rlist;
  vector<OERing*>::iterator i;

  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();
  if (!mol->HasSSSRPerceived())
    mol->FindSSSR();

 if (!((OEAtom*)this)->HasFlag(OE_RING_ATOM)) return(false);

  rlist = mol->GetSSSR();
  for (i = rlist.begin();i != rlist.end();i++)
    if ((*i)->IsInRing(GetIdx()) && (*i)->PathSize() == size)
      return(true);
  
  return(false);
}

unsigned int OEAtom::MemberOfRingCount() const
{
  vector<OERing*> rlist;
  vector<OERing*>::iterator i;

  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();

  if (!mol->HasSSSRPerceived())
    mol->FindSSSR();

  if (!((OEAtom*)this)->IsInRing()) return(0);

  int count=0;
  rlist = mol->GetSSSR();

  for (i = rlist.begin();i != rlist.end();i++)
    if ((*i)->IsInRing(GetIdx()))
      count++;
  
  return((unsigned int)count);
}

unsigned int OEAtom::CountFreeOxygens() const
{
  unsigned int count = 0;
  OEAtom *atom;
  vector<OEBond*>::iterator i;

  for (atom = ((OEAtom*)this)->BeginNbrAtom(i);atom;atom = ((OEAtom*)this)->NextNbrAtom(i))
	  if (atom->IsOxygen() && atom->GetHvyValence() == 1)
		  count++;

  return(count);
}

unsigned int OEAtom::BOSum() const
{
  unsigned int bo;
  unsigned int bosum=0;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  
  for (bond = ((OEAtom*)this)->BeginBond(i);bond;bond = ((OEAtom*)this)->NextBond(i))
    {
      bo = bond->GetBO();
      bosum += (bo < 4) ? 2*bo : 3;
    }

  bosum /= 2;
  return(bosum);
}

unsigned int OEAtom::KBOSum() const
{
  OEBond *bond;
  unsigned int bosum = 0;
  vector<OEBond*>::iterator i;
  
  bosum = GetImplicitValence();

  for (bond = ((OEAtom*)this)->BeginBond(i);bond;bond = ((OEAtom*)this)->NextBond(i))
    {
      if (bond->IsKDouble()) bosum++;
      else if (bond->IsKTriple()) bosum += 2;
    }

  return(bosum);
}

unsigned int OEAtom::ImplicitHydrogenCount() const
     //handles H,C,N,S,O,X
{
  OEMol *mol = (OEMol*)((OEAtom*)this)->GetParent();
  if (mol && !mol->HasImplicitValencePerceived())
    atomtyper.AssignImplicitValence(*((OEMol*)((OEAtom*)this)->GetParent()));

  int impval = _impval - GetHvyValence();
  return((impval>0)?impval:0);
}

unsigned int OEAtom::ExplicitHydrogenCount() const
{
  int numH=0;
  OEAtom *atom;
  vector<OEBond*>::iterator i;
  for (atom = ((OEAtom*)this)->BeginNbrAtom(i);atom;atom = ((OEAtom*)this)->NextNbrAtom(i))
    if (atom->IsHydrogen())
      numH++;

  return(numH);
}

bool OEAtom::DeleteBond(OEBond *bond)
{
  vector<OEEdgeBase*>::iterator i;
  for (i = _vbond.begin();i != _vbond.end();i++)
    if ((OEBond*)bond == *i)
    {
      _vbond.erase(i);
      return(true);
    }
  return(false);
}

OEBond *OEAtom::BeginBond(vector<OEBond*>::iterator &i) 
{
	i = (vector<OEBond*>::iterator)_vbond.begin();
	return((i == (vector<OEBond*>::iterator)_vbond.end()) ? (OEBond*)NULL : (OEBond*)*i);
}

OEBond *OEAtom::NextBond(vector<OEBond*>::iterator &i) 
{
	i++;
	return((i == (vector<OEBond*>::iterator)_vbond.end()) ? (OEBond*)NULL : (OEBond*)*i);
}

OEAtom *OEAtom::BeginNbrAtom(vector<OEBond*>::iterator &i)
{
	i = (vector<OEBond*>::iterator)_vbond.begin();
  return((i != (vector<OEBond*>::iterator)_vbond.end()) ? (*i)->GetNbrAtom(this):NULL);
}

OEAtom *OEAtom::NextNbrAtom(vector<OEBond*>::iterator &i) 
{
  i++; 
  return((i != (vector<OEBond*>::iterator)_vbond.end()) ? (*i)->GetNbrAtom(this):NULL);
}


bool OEAtom::GetNewBondVector(Vector &v,float length)
{
  // ***experimental code***

  OEAtom *atom;
  vector<OEBond*>::iterator i,j;
  v = VZero;
  if (GetValence() == 0)
    {
      v = VX;
      v *= length;
      v += GetVector();
      return(true);
    }

  if (GetValence() == 1)
    {
      Vector vtmp,v1,v2;
      atom = BeginNbrAtom(i);
      if (atom)
	vtmp = GetVector() - atom->GetVector();

      if (GetHyb() == 2 || (IsOxygen() && HasAlphaBetaUnsat()))
	{
	  bool quit = false;
	  OEAtom *a1,*a2;
	  v2 = VZero;
	  for (a1 = BeginNbrAtom(i);a1 && !quit;a1 = NextNbrAtom(i))
	    for (a2 = a1->BeginNbrAtom(j);a2 && !quit;a2 = a1->NextNbrAtom(j))
	      if (a1 && a2 && a2 != this)
		{
		  v2 = a1->GetVector() - a2->GetVector();
		  quit = true;
		}

	  if (v2 == VZero)
	    {
	      v1 = cross(vtmp,VX);
	      v2 = cross(vtmp,VY);
	      if (v1.length() < v2.length()) v1 = v2;
	    }
	  else
	      v1 = cross(vtmp,v2);

	  Matrix3x3 m;
	  m.RotAboutAxisByAngle(v1,60.0);
	  v = m*vtmp;
	  v.normalize();
	}

      if (GetHyb() == 3)
	{
	  v1 = cross(vtmp,VX);
	  v2 = cross(vtmp,VY);
	  if (v1.length() < v2.length()) v1 = v2;
	  Matrix3x3 m;
	  m.RotAboutAxisByAngle(v1,70.5);
	  v = m*vtmp;
	  v.normalize();
	}
      if (GetHyb() == 1) v = vtmp;

      v *= length;
      v += GetVector();
      return(true);
    }

  if (GetValence() == 2)
    {
      Vector v1,v2,vtmp,vsum,vnorm;
      atom = BeginNbrAtom(i);if (!atom) return(false);
      v1 = GetVector() - atom->GetVector();
      atom = NextNbrAtom(i); if (!atom) return(false);
      v2 = GetVector() - atom->GetVector();
      v1.normalize();v2.normalize();
      vsum = v1+v2;
      vsum.normalize();

      if (GetHyb() == 2) v = vsum;

      if (GetHyb() == 3)
	{
	  vnorm = cross(v2,v1);
	  vnorm.normalize();

#ifndef ONE_OVER_SQRT3
#define ONE_OVER_SQRT3  0.577350269f
#endif //SQRT_TWO_THIRDS
#ifndef SQRT_TWO_THIRDS
#define SQRT_TWO_THIRDS 0.816496581f
#endif //ONE_OVER_SQRT3

	  vsum *= ONE_OVER_SQRT3;
	  vnorm *= SQRT_TWO_THIRDS;

#undef ONE_OVER_SQRT3
#undef SQRT_TWO_THIRDS

	  v = vsum + vnorm;
	}
      v *= length;
      v += GetVector();
      return(true);
    }

  if (GetValence() == 3)
    {
      Vector vtmp,vsum;
      OEAtom *atom;
      vector<OEBond*>::iterator i;
      for (atom = BeginNbrAtom(i);atom;atom = NextNbrAtom(i))
	{
	  vtmp = GetVector() - atom->GetVector();
	  vtmp.normalize();
	  vtmp /= 3.0;
	  vsum += vtmp;
	}
      vsum.normalize();
      v = vsum;
      v *= length;
      v += GetVector();
      return(true);
    }

  return(true);
}

bool OEAtom::HtoMethyl()
{
  if (!IsHydrogen()) return(false);
  OEMol *mol = (OEMol*)GetParent();

  mol->BeginModify();

  SetAtomicNum(6); SetType("C3"); SetHyb(3);

  OEAtom *atom;
  OEBond *bond;
  vector<OEBond*>::iterator i;
  atom  = BeginNbrAtom(i);
  bond = *i;
  if (!atom)
    {
      mol->EndModify();
      return(false);
    }

  float br1,br2;
  br1 = etab.CorrectedBondRad(6,3);
  br2 = etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
  bond->SetLength(atom,br1+br2);

  OEAtom *hatom;
  br2 = etab.CorrectedBondRad(1,0);
  Vector v;
  
  for (int j = 0;j < 3;j++)
    {
      hatom = mol->NewAtom();
      hatom->SetAtomicNum(1);
      hatom->SetType("H");

      GetNewBondVector(v,br1+br2); 
      hatom->SetVector(v);
      mol->AddBond(GetIdx(),mol->NumAtoms(),1);
    }

  mol->EndModify();
  return(true);
}

static void ApplyRotMatToBond(OEMol &mol,Matrix3x3 &m,OEAtom *a1,OEAtom *a2)
{
  vector<int> children;
  mol.FindChildren(children,a1->GetIdx(),a2->GetIdx());
  children.push_back(a2->GetIdx());

  Vector v;
  vector<int>::iterator i;
  for (i = children.begin();i != children.end();i++)
    {
      v = mol.GetAtom(*i)->GetVector();
      v -= a1->GetVector();
      v *= m;
      v += a1->GetVector();
      mol.GetAtom(*i)->SetVector(v);
    }

}

bool OEAtom::SetHybAndGeom(int hyb)
{
  //if (hyb == GetHyb()) return(true);
  if (hyb == 0 && GetHvyValence() > 1) return(false);
  if (hyb == 1 && GetHvyValence() > 2) return(false);
  if (hyb == 2 && GetHvyValence() > 3) return(false);
  if (hyb == 3 && GetHvyValence() > 4) return(false);

  OEMol *mol = (OEMol*)GetParent();

  OEAtom *nbr;
  vector<OEAtom*> delatm;
  vector<OEBond*>::iterator i;

  for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
    if (nbr->IsHydrogen())
      delatm.push_back(nbr);
  
  //delete attached hydrogens
  mol->IncrementMod();
  vector<OEAtom*>::iterator j;
  for (j = delatm.begin();j != delatm.end();j++) mol->DeleteAtom(*j);
  mol->DecrementMod();

  float targetAngle;
  if (hyb == 3)     targetAngle = 109.5;
  else if (hyb == 2) targetAngle = 120.0;
  else if (hyb == 1) targetAngle = 180.0;
  else               targetAngle = 0.0;

  //adjust attached acyclic bond lengths
  float br1,br2;
  br1 = etab.CorrectedBondRad(GetAtomicNum(),hyb);
  for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
    if (!(*i)->IsInRing())
    {
      br2 = etab.CorrectedBondRad(nbr->GetAtomicNum(),nbr->GetHyb());
      (*i)->SetLength(this,br1+br2);
    }

  if (GetValence() > 1)
    {
      float angle;
      Matrix3x3 m;
      Vector v1,v2,v3,v4,n,s;
      OEAtom *r1,*r2,*r3,*a1,*a2,*a3,*a4;
      r1 = r2 = r3 = a1 = a2 = a3 = a4 = NULL;

      //find ring atoms first
      for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
	if ((*i)->IsInRing())
	  if (!r1)     r1 = nbr;
	  else if (!r2) r2 = nbr;
	  else if (!r3) r3 = nbr;

      //find non-ring atoms
      for (nbr = BeginNbrAtom(i);nbr;nbr = NextNbrAtom(i))
	if (!(*i)->IsInRing())
	  if (!a1)     a1 = nbr;
	  else if (!a2) a2 = nbr;
	  else if (!a3) a3 = nbr;
	  else if (!a4) a4 = nbr;

      //adjust geometries of heavy atoms according to hybridization
      if (hyb == 1)
	{
	  if (a2)
	    {
	      v1 = a1->GetVector()-GetVector();v1.normalize();
	      v2 = a2->GetVector()-GetVector();v2.normalize();
	      n = cross(v1,v2);
	      angle = VectorAngle(v1,v2)-targetAngle;
	      m.RotAboutAxisByAngle(n,-angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	}
      else if (hyb == 2)
	{
	  if (r1 && r2 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = a1->GetVector()-GetVector();
	      s = v1+v2; s.normalize(); s *= -1.0f;
	      n = cross(s,v3);
	      angle = VectorAngle(s,v3);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	  else
	    {
	      if (a2)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  n = cross(v1,v2);
		  angle = VectorAngle(v1,v2)-targetAngle;
		  m.RotAboutAxisByAngle(n,-angle);
		  ApplyRotMatToBond(*mol,m,this,a1);
		}
	      if (a3)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  v3 = a3->GetVector()-GetVector();
		  s = v1+v2; s.normalize(); s *= -1.0f;
		  n = cross(s,v3);
		  angle = VectorAngle(s,v3);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a3);
		}
	    }
	}
      else if (hyb == 3)
	{
	  if (r1 && r2 && r3 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = r3->GetVector()-GetVector();v3.normalize();
	      v4 = a1->GetVector()-GetVector();
	      s = v1 + v2 + v3; s *= -1.0f; s.normalize();
	      n = cross(s,v4);
	      angle = VectorAngle(s,v4);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);
	    }
	  else if (r1 && r2 && !r3 && a1)
	    {
	      v1 = r1->GetVector()-GetVector();v1.normalize();
	      v2 = r2->GetVector()-GetVector();v2.normalize();
	      v3 = a1->GetVector()-GetVector();
	      s = v1+v2; s *= -1.0f; s.normalize(); 
	      n = cross(v1,v2);      n.normalize();
	      s *= 0.577350269f; //1/sqrt(3)
	      n *= 0.816496581f; //sqrt(2/3)
	      s += n; s.normalize();
	      n = cross(s,v3);
	      angle = VectorAngle(s,v3);
	      m.RotAboutAxisByAngle(n,angle);
	      ApplyRotMatToBond(*mol,m,this,a1);

	      if (a2)
		{
		  v1 = r1->GetVector()-GetVector();v1.normalize();
		  v2 = r2->GetVector()-GetVector();v2.normalize();
		  v3 = a1->GetVector()-GetVector();v3.normalize();
		  v4 = a2->GetVector()-GetVector();
		  s = v1 + v2 + v3; s *= -1.0f; s.normalize();
		  n = cross(s,v4);
		  angle = VectorAngle(s,v4);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a2);
		}
	    }
	  else
	    {
	      if (a2)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  n = cross(v1,v2);
		  angle = VectorAngle(v1,v2)-targetAngle;
		  m.RotAboutAxisByAngle(n,-angle);
		  ApplyRotMatToBond(*mol,m,this,a1);
		}
	      if (a3)
		{
		  v1 = a1->GetVector()-GetVector();v1.normalize();
		  v2 = a2->GetVector()-GetVector();v2.normalize();
		  v3 = a3->GetVector()-GetVector();
		  s = v1+v2; s *= -1.0f; s.normalize(); 
		  n = cross(v1,v2);      n.normalize();
		  s *= 0.577350269f; //1/sqrt(3)
		  n *= 0.816496581f; //sqrt(2/3)
		  s += n; s.normalize();
		  n = cross(s,v3);
		  angle = VectorAngle(s,v3);
		  m.RotAboutAxisByAngle(n,angle);
		  ApplyRotMatToBond(*mol,m,this,a3);
		}
	    }
	}
    }

  //add hydrogens back to atom
  int impval=1;
  switch (GetAtomicNum())
    {
    case 6: 
      if (hyb == 3) impval = 4;
      if (hyb == 2) impval = 3;
      if (hyb == 1) impval = 2;
      break;
    case 7: 
      if (hyb == 3) impval = 3;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 1;
      break;
    case 8: 
      if (hyb == 3) impval = 2;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 2;
    case 16: 
      if (hyb == 3) impval = 2;
      if (hyb == 2) impval = 2;
      if (hyb == 1) impval = 0;
      break;
    case 15: 
      if (hyb == 3) impval = 4;
      if (hyb == 2) impval = 3;
      if (hyb == 1) impval = 2;
      break;
    default:
      impval = 1;
    }

  int hcount = impval-GetHvyValence();
  if (hcount)
    {
      int k;
      Vector v;
      OEAtom *atom;
      float brsum = etab.CorrectedBondRad(1,0)+
	etab.CorrectedBondRad(GetAtomicNum(),GetHyb());
      SetHyb(hyb);

      mol->BeginModify();
      for (k = 0;k < hcount;k++)
	{
	  GetNewBondVector(v,brsum);
	  atom = mol->NewAtom();
	  atom->SetAtomicNum(1);
	  atom->SetType("H");
	  atom->SetVector(v);
	  mol->AddBond(atom->GetIdx(),GetIdx(),1);
	}
      mol->EndModify();
    }


  return(true);
}

OEBond *OEAtom::GetBond(OEAtom *nbr)
{
  OEBond *bond;
  vector<OEBond *>::iterator i;
  for (bond = BeginBond(i) ; bond ; bond = NextBond(i))
    if (bond->GetNbrAtom(this) == nbr)
      return bond;
  return NULL;
}

}
