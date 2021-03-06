#include "obtest.h"
#include <openbabel/tautomer.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace OpenBabel;




void testEnumerateTautomers(const std::string &smiles, int numTautomers)
{
  class Functor : public TautomerFunctor
  {
    public:
      int numTautomers;

      Functor() : numTautomers(0) {}
      void operator()(OBMol*)
      {
        numTautomers++;
      }
  };

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.ReadString(&mol, smiles);

  Functor functor;
  EnumerateTautomers(&mol, functor);

  OB_COMPARE( functor.numTautomers, numTautomers );
}

void testCanonicalTautomers(const std::string &smiles)
{
  class Functor : public TautomerFunctor
  {
    public:
      std::vector<std::string> tautomers;

      void operator()(OBMol *mol)
      {
        OBConversion conv;
        conv.SetOutFormat("can");
        tautomers.push_back(conv.WriteString(mol));
      }
  };

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");
  conv.ReadString(&mol, smiles);

  // enumerate all tautomers
  Functor functor;
  EnumerateTautomers(&mol, functor);

  // check to make sure all tautomers result in the same canonical tautomer
  const std::vector<std::string> &tautomers = functor.tautomers;
  if (tautomers.empty())
    return;

  std::vector<std::string> canonicalTautomers;
  for (std::size_t i = 0; i < tautomers.size(); ++i) {
    OBMol mol2;
    conv.ReadString(&mol2, tautomers[i]);
    CanonicalTautomer(&mol2);
    canonicalTautomers.push_back(conv.WriteString(&mol2));
  }

  canonicalTautomers.erase(std::unique(canonicalTautomers.begin(), canonicalTautomers.end()), canonicalTautomers.end());

  OB_COMPARE( canonicalTautomers.size(), 1 );
}




int main()
{
  // guanine
  testEnumerateTautomers("Nc1nc2ncnc2c([nH]1)O", 15);
  testCanonicalTautomers("Nc1nc2ncnc2c([nH]1)O");


}
