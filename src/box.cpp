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

using namespace std;

namespace OpenBabel 
{

bool ReadBox(istream &ifs, OBMol &mol,const char *title)
{
  char buffer[BUFF_SIZE];
  vector<string> vs;
  vector<string>::iterator i;
  OBAtom atom;

  mol.BeginModify();
  
  while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"END",3))
    {
      if (EQn(buffer,"ATOM",4))
	{
	  string sbuf = &buffer[6]; 
	  /* X, Y, Z */
	  string x = sbuf.substr(24,8);
	  string y = sbuf.substr(32,8);
	  string z = sbuf.substr(40,8);
	  Vector v(atof(x.c_str()),atof(y.c_str()),atof(z.c_str()));
	  atom.SetVector(v);
	  if (!mol.AddAtom(atom)) return(false);
	}
    
      if (EQn(buffer,"CONECT",6))
	{
	  tokenize(vs,buffer);
	  if (!vs.empty() && vs.size() > 2)
	    for (i = vs.begin(),i+=2;i != vs.end();i++)
	      mol.AddBond(atoi(vs[1].c_str()),atoi((*i).c_str()),1);
	}
    }

  mol.EndModify();
  mol.SetTitle(title);
  return(true);
}

bool WriteBox(ostream &ofs,OBMol &mol,float margin)
{
  char buffer[BUFF_SIZE];
  Vector vcenter,vmin,vmax,vmid,vdim;
  
  OBAtom *atom;
  vector<OBNodeBase*>::iterator i;
  vmax.Set(-10E10,-10E10,-10E10);
  vmin.Set( 10E10, 10E10, 10E10);

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      vcenter += atom->GetVector();
      if (atom->x() < vmin.x()) vmin.SetX(atom->x());
      if (atom->y() < vmin.y()) vmin.SetY(atom->y());
      if (atom->z() < vmin.z()) vmin.SetZ(atom->z());

      if (atom->x() > vmax.x()) vmax.SetX(atom->x());
      if (atom->y() > vmax.y()) vmax.SetY(atom->y());
      if (atom->z() > vmax.z()) vmax.SetZ(atom->z());
    }
  vcenter /= (float)mol.NumAtoms();

  Vector vmarg(margin,margin,margin);
  vmin -= vmarg; vmax += vmarg;
  vdim = vmax - vmin;
  vmid = vmin+vmax; vmid /= 2.0f;

  ofs << "HEADER    CORNERS OF BOX" << endl;
  sprintf(buffer,"REMARK    CENTER (X Y Z)      %10.3f %10.3f %10.3f",
	  vmid.x(),vmid.y(),vmid.z()); ofs << buffer << endl;
  sprintf(buffer,"REMARK    DIMENSIONS (X Y Z)  %10.3f %10.3f %10.3f",
	  vdim.x(),vdim.y(),vdim.z()); ofs << buffer << endl;
  vdim /= 2.0f;

  Vector vtmp;
  int j;
  for (j = 1;j <= 8;j++)
    {
      switch(j)
	{
	case 1: vtmp = vmid-vdim; break;
	case 2: vtmp.SetX(vmid.x()+vdim.x()); break;
	case 3: vtmp.SetZ(vmid.z()+vdim.z()); break;
	case 4: vtmp.SetX(vmid.x()-vdim.x()); break;
	case 5: 
	  vtmp = vmid-vdim; vtmp.SetY(vmid.y()+vdim.y());
	  break;
	case 6: 
	  vtmp = vmid+vdim; vtmp.SetZ(vmid.z()-vdim.z());
	  break;
	case 7: vtmp = vmid+vdim; break;
	case 8: vtmp.SetX(vmid.x()-vdim.x()); break;
	}
      sprintf(buffer,"ATOM      %d  DUA BOX     1    %8.3f%8.3f%8.3f",
	      j,vtmp.x(),vtmp.y(),vtmp.z());
      ofs << buffer << endl;
    }

  ofs << "CONECT    1    2    4    5" << endl;
  ofs << "CONECT    2    1    3    6" << endl;
  ofs << "CONECT    3    2    4    7" << endl;
  ofs << "CONECT    4    1    3    8" << endl;
  ofs << "CONECT    5    1    6    8" << endl;
  ofs << "CONECT    6    2    5    7" << endl;
  ofs << "CONECT    7    3    6    8" << endl;
  ofs << "CONECT    8    4    5    7" << endl;
  
  return(true);
}

}
