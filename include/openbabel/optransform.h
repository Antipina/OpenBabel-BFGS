/**********************************************************************
optransform.h: makes option to transform molecule as specified in a datafile
Copyright (C) 2008 Chris Morley
 
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
#include <openbabel/mol.h>
#include <openbabel/phmodel.h>
#include <openbabel/op.h>
#include <vector>

namespace OpenBabel
{
  ///Option that applies molecular transforms (OBChemTsfm class) read from a datafile
class OpTransform : public OBOp
{
public:
  //! constructor. Each instance provides an ID, a datafile and a description.
  OpTransform(const char* ID, const char* filename, const char* descr)
    : OBOp(ID, false), _filename(filename), _descr(descr), _dataLoaded(false){}

  ~OpTransform(){}

  virtual const char* Description();

  //!Checks that this op is being applied to the right kind of object(OBMol)
  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }

  //! Carries out the transform
  virtual bool Do(OBBase* pOb, OpMap* pOptions=NULL, const char* OptionText=NULL);

  virtual OpTransform* MakeInstance(const std::vector<std::string>& textlines)
{
  OpTransform* pTransform = new OpTransform(
    textlines[1].c_str(),textlines[2].c_str(),textlines[3].c_str());
  pTransform->_textlines = textlines;
  return pTransform;
}

private:
  bool Initialize();
  void ParseLine(const char *buffer);

private:
  const char* _filename;
  const char* _descr;
  std::vector<std::string> _textlines;

  bool _dataLoaded;
  std::vector<OBChemTsfm> _transforms;
};

}//namespace

