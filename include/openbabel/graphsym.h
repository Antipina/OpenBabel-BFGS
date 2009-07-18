/**********************************************************************
graphsym.h - Class for handling graph symmetry.

To determine copyright, please analyse the Subversion commit log.

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

#ifndef OB_GRAPHSYM_H
#define OB_GRAPHSYM_H

#include <openbabel/babelconfig.h>

#ifndef EXTERN
#  define EXTERN extern
#endif

using namespace std;

namespace OpenBabel {

  class OBBitVec;
  class OBMol;
  class OBAtom;
  class OBMol;

  class OBAPI OBGraphSym {
    private:
      OBBitVec* _frag_atoms;
      OBMol* _pmol;
      std::vector<int> _gtd;

      static bool CompareUnsigned(const unsigned int &a,const unsigned int &b);
      static bool ComparePairFirst(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b);
      static bool ComparePairSecond(const std::pair<OBAtom*,unsigned int> &a,const std::pair<OBAtom*,unsigned int> &b);
      static bool CompareBondPairSecond(const std::pair<OBBond*,unsigned int> &a,const std::pair<OBBond*,unsigned int> &b);

      unsigned int GetValence(OBAtom *atom);
      unsigned int GetHvyValence(OBAtom *atom);
      unsigned int GetHvyBondSum(OBAtom *atom);

      bool CalcGTDVector();
      void FindRingAtoms(OBBitVec &ring_atoms);
      void CreateNewClassVector(std::vector<std::pair<OBAtom*,unsigned int> > &vp1,
                                std::vector<std::pair<OBAtom*,unsigned int> > &vp2,
                                int natoms);
      void GetGIVector(std::vector<unsigned int> &vid);
      void CountAndRenumberClasses(std::vector<std::pair<OBAtom*,unsigned int> > &vp, unsigned int &count);
      int ExtendInvariants(std::vector<std::pair<OBAtom*, unsigned int> > &symmetry_classes,
                           int nfragatoms, int natoms);
      int CalculateSymmetry(std::vector<std::pair<OBAtom*, unsigned int> > &symmetry_classes);

    public:
      //! Constructor
      OBGraphSym(OBMol* pmol = NULL, OBBitVec* frag_atoms = NULL)
      {
        _pmol = pmol;
        if (_pmol) {
          if (frag_atoms)
            _frag_atoms = frag_atoms;
          else
          {
            _frag_atoms = new OBBitVec(_pmol->NumAtoms());
            FOR_ATOMS_OF_MOL(a, _pmol)
              _frag_atoms->SetBitOn(a->GetIdx());
          }
        }
        else // No molecule supplied
          _frag_atoms = NULL;
      }
      //OBGraphSym(const OBMol&);
      //OBGraphSym(const OBMol&, const OBBitVec&);
      //! Destructor
      virtual ~OBGraphSym();

      bool GetGTDVector(std::vector<int> &gtd) {
        bool success = CalcGTDVector();
        if (success)
          gtd =_gtd; // Copy _gtd into gtd
        return success;
      }
    };

      //&OBBitVec GetFragment();
      //SetFragment(&OBBitVec);
} // namespace OpenBabel

//! \file graphsym.h
//! \brief XXXX

  #endif // OB_GRAPHSYM_H