/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMol extends OBGraphBase {
  private long swigCPtr;

  protected OBMol(long cPtr, boolean cMemoryOwn) {
    super(net.sourceforge.openbabelJNI.SWIGOBMolUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMol obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_OBMol(swigCPtr);
    }
    swigCPtr = 0;
    super.delete();
  }

  public OBMol() {
    this(net.sourceforge.openbabelJNI.new_OBMol__SWIG_0(), true);
  }

  public OBMol(OBMol arg0) {
    this(net.sourceforge.openbabelJNI.new_OBMol__SWIG_1(OBMol.getCPtr(arg0), arg0), true);
  }

  public void ReserveAtoms(int natoms) {
    net.sourceforge.openbabelJNI.OBMol_ReserveAtoms(swigCPtr, this, natoms);
  }

  public OBAtom CreateAtom() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_CreateAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBBond CreateBond() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_CreateBond(swigCPtr, this);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBResidue CreateResidue() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_CreateResidue(swigCPtr, this);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public void DestroyAtom(OBNodeBase arg0) {
    net.sourceforge.openbabelJNI.OBMol_DestroyAtom(swigCPtr, this, OBNodeBase.getCPtr(arg0), arg0);
  }

  public void DestroyBond(OBEdgeBase arg0) {
    net.sourceforge.openbabelJNI.OBMol_DestroyBond(swigCPtr, this, OBEdgeBase.getCPtr(arg0), arg0);
  }

  public void DestroyResidue(OBResidue arg0) {
    net.sourceforge.openbabelJNI.OBMol_DestroyResidue(swigCPtr, this, OBResidue.getCPtr(arg0), arg0);
  }

  public boolean AddAtom(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_AddAtom(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean AddBond(int arg0, int arg1, int arg2, int flags, int insertpos) {
    return net.sourceforge.openbabelJNI.OBMol_AddBond__SWIG_0(swigCPtr, this, arg0, arg1, arg2, flags, insertpos);
  }

  public boolean AddBond(int arg0, int arg1, int arg2, int flags) {
    return net.sourceforge.openbabelJNI.OBMol_AddBond__SWIG_1(swigCPtr, this, arg0, arg1, arg2, flags);
  }

  public boolean AddBond(int arg0, int arg1, int arg2) {
    return net.sourceforge.openbabelJNI.OBMol_AddBond__SWIG_2(swigCPtr, this, arg0, arg1, arg2);
  }

  public boolean AddBond(OBBond arg0) {
    return net.sourceforge.openbabelJNI.OBMol_AddBond__SWIG_3(swigCPtr, this, OBBond.getCPtr(arg0), arg0);
  }

  public boolean AddResidue(OBResidue arg0) {
    return net.sourceforge.openbabelJNI.OBMol_AddResidue(swigCPtr, this, OBResidue.getCPtr(arg0), arg0);
  }

  public boolean InsertAtom(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_InsertAtom(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean DeleteAtom(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_DeleteAtom(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean DeleteBond(OBBond arg0) {
    return net.sourceforge.openbabelJNI.OBMol_DeleteBond(swigCPtr, this, OBBond.getCPtr(arg0), arg0);
  }

  public boolean DeleteResidue(OBResidue arg0) {
    return net.sourceforge.openbabelJNI.OBMol_DeleteResidue(swigCPtr, this, OBResidue.getCPtr(arg0), arg0);
  }

  public OBAtom NewAtom() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NewAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBResidue NewResidue() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NewResidue(swigCPtr, this);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public void BeginModify() {
    net.sourceforge.openbabelJNI.OBMol_BeginModify(swigCPtr, this);
  }

  public void EndModify(boolean nukePerceivedData) {
    net.sourceforge.openbabelJNI.OBMol_EndModify__SWIG_0(swigCPtr, this, nukePerceivedData);
  }

  public void EndModify() {
    net.sourceforge.openbabelJNI.OBMol_EndModify__SWIG_1(swigCPtr, this);
  }

  public int GetMod() {
    return net.sourceforge.openbabelJNI.OBMol_GetMod(swigCPtr, this);
  }

  public void IncrementMod() {
    net.sourceforge.openbabelJNI.OBMol_IncrementMod(swigCPtr, this);
  }

  public void DecrementMod() {
    net.sourceforge.openbabelJNI.OBMol_DecrementMod(swigCPtr, this);
  }

  public int GetFlags() {
    return net.sourceforge.openbabelJNI.OBMol_GetFlags(swigCPtr, this);
  }

  public String GetTitle() {
    return net.sourceforge.openbabelJNI.OBMol_GetTitle(swigCPtr, this);
  }

  public long NumAtoms() {
    return net.sourceforge.openbabelJNI.OBMol_NumAtoms(swigCPtr, this);
  }

  public long NumBonds() {
    return net.sourceforge.openbabelJNI.OBMol_NumBonds(swigCPtr, this);
  }

  public long NumHvyAtoms() {
    return net.sourceforge.openbabelJNI.OBMol_NumHvyAtoms(swigCPtr, this);
  }

  public long NumResidues() {
    return net.sourceforge.openbabelJNI.OBMol_NumResidues(swigCPtr, this);
  }

  public long NumRotors() {
    return net.sourceforge.openbabelJNI.OBMol_NumRotors(swigCPtr, this);
  }

  public OBAtom GetAtom(int arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetAtom(swigCPtr, this, arg0);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom GetFirstAtom() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetFirstAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBBond GetBond(int arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetBond__SWIG_0(swigCPtr, this, arg0);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond GetBond(int arg0, int arg1) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetBond__SWIG_1(swigCPtr, this, arg0, arg1);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond GetBond(OBAtom bgn, OBAtom end) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetBond__SWIG_2(swigCPtr, this, OBAtom.getCPtr(bgn), bgn, OBAtom.getCPtr(end), end);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBResidue GetResidue(int arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetResidue(swigCPtr, this, arg0);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t GetInternalCoord() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t(net.sourceforge.openbabelJNI.OBMol_GetInternalCoord(swigCPtr, this), true);
  }

  public double GetTorsion(int arg0, int arg1, int arg2, int arg3) {
    return net.sourceforge.openbabelJNI.OBMol_GetTorsion__SWIG_0(swigCPtr, this, arg0, arg1, arg2, arg3);
  }

  public double GetTorsion(OBAtom a, OBAtom b, OBAtom c, OBAtom d) {
    return net.sourceforge.openbabelJNI.OBMol_GetTorsion__SWIG_1(swigCPtr, this, OBAtom.getCPtr(a), a, OBAtom.getCPtr(b), b, OBAtom.getCPtr(c), c, OBAtom.getCPtr(d), d);
  }

  public double GetAngle(OBAtom a, OBAtom b, OBAtom c) {
    return net.sourceforge.openbabelJNI.OBMol_GetAngle(swigCPtr, this, OBAtom.getCPtr(a), a, OBAtom.getCPtr(b), b, OBAtom.getCPtr(c), c);
  }

  public String GetFormula() {
    return net.sourceforge.openbabelJNI.OBMol_GetFormula(swigCPtr, this);
  }

  public String GetSpacedFormula(int ones, String sp) {
    return net.sourceforge.openbabelJNI.OBMol_GetSpacedFormula__SWIG_0(swigCPtr, this, ones, sp);
  }

  public String GetSpacedFormula(int ones) {
    return net.sourceforge.openbabelJNI.OBMol_GetSpacedFormula__SWIG_1(swigCPtr, this, ones);
  }

  public String GetSpacedFormula() {
    return net.sourceforge.openbabelJNI.OBMol_GetSpacedFormula__SWIG_2(swigCPtr, this);
  }

  public double GetEnergy() {
    return net.sourceforge.openbabelJNI.OBMol_GetEnergy(swigCPtr, this);
  }

  public double GetMolWt() {
    return net.sourceforge.openbabelJNI.OBMol_GetMolWt(swigCPtr, this);
  }

  public double GetExactMass() {
    return net.sourceforge.openbabelJNI.OBMol_GetExactMass(swigCPtr, this);
  }

  public int GetTotalCharge() {
    return net.sourceforge.openbabelJNI.OBMol_GetTotalCharge(swigCPtr, this);
  }

  public long GetTotalSpinMultiplicity() {
    return net.sourceforge.openbabelJNI.OBMol_GetTotalSpinMultiplicity(swigCPtr, this);
  }

  public int GetDimension() {
    return net.sourceforge.openbabelJNI.OBMol_GetDimension(swigCPtr, this);
  }

  public SWIGTYPE_p_double GetCoordinates() {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetCoordinates(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBRing_p_t GetSSSR() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBRing_p_t(net.sourceforge.openbabelJNI.OBMol_GetSSSR(swigCPtr, this), false);
  }

  public boolean AutomaticFormalCharge() {
    return net.sourceforge.openbabelJNI.OBMol_AutomaticFormalCharge(swigCPtr, this);
  }

  public boolean AutomaticPartialCharge() {
    return net.sourceforge.openbabelJNI.OBMol_AutomaticPartialCharge(swigCPtr, this);
  }

  public void SetTitle(String title) {
    net.sourceforge.openbabelJNI.OBMol_SetTitle__SWIG_0(swigCPtr, this, title);
  }

  public void SetTitle(SWIGTYPE_p_std__string title) {
    net.sourceforge.openbabelJNI.OBMol_SetTitle__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__string.getCPtr(title));
  }

  public void SetFormula(String molFormula) {
    net.sourceforge.openbabelJNI.OBMol_SetFormula(swigCPtr, this, molFormula);
  }

  public void SetEnergy(double energy) {
    net.sourceforge.openbabelJNI.OBMol_SetEnergy(swigCPtr, this, energy);
  }

  public void SetDimension(int d) {
    net.sourceforge.openbabelJNI.OBMol_SetDimension(swigCPtr, this, d);
  }

  public void SetTotalCharge(int charge) {
    net.sourceforge.openbabelJNI.OBMol_SetTotalCharge(swigCPtr, this, charge);
  }

  public void SetTotalSpinMultiplicity(long spin) {
    net.sourceforge.openbabelJNI.OBMol_SetTotalSpinMultiplicity(swigCPtr, this, spin);
  }

  public void SetInternalCoord(SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t int_coord) {
    net.sourceforge.openbabelJNI.OBMol_SetInternalCoord(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t.getCPtr(int_coord));
  }

  public void SetAutomaticFormalCharge(boolean val) {
    net.sourceforge.openbabelJNI.OBMol_SetAutomaticFormalCharge(swigCPtr, this, val);
  }

  public void SetAutomaticPartialCharge(boolean val) {
    net.sourceforge.openbabelJNI.OBMol_SetAutomaticPartialCharge(swigCPtr, this, val);
  }

  public void SetAromaticPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetAromaticPerceived(swigCPtr, this);
  }

  public void SetSSSRPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetSSSRPerceived(swigCPtr, this);
  }

  public void SetRingAtomsAndBondsPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetRingAtomsAndBondsPerceived(swigCPtr, this);
  }

  public void SetAtomTypesPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetAtomTypesPerceived(swigCPtr, this);
  }

  public void SetChainsPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetChainsPerceived(swigCPtr, this);
  }

  public void SetChiralityPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetChiralityPerceived(swigCPtr, this);
  }

  public void SetPartialChargesPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetPartialChargesPerceived(swigCPtr, this);
  }

  public void SetHybridizationPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetHybridizationPerceived(swigCPtr, this);
  }

  public void SetImplicitValencePerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetImplicitValencePerceived(swigCPtr, this);
  }

  public void SetKekulePerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetKekulePerceived(swigCPtr, this);
  }

  public void SetClosureBondsPerceived() {
    net.sourceforge.openbabelJNI.OBMol_SetClosureBondsPerceived(swigCPtr, this);
  }

  public void SetHydrogensAdded() {
    net.sourceforge.openbabelJNI.OBMol_SetHydrogensAdded(swigCPtr, this);
  }

  public void SetCorrectedForPH() {
    net.sourceforge.openbabelJNI.OBMol_SetCorrectedForPH(swigCPtr, this);
  }

  public void SetAromaticCorrected() {
    net.sourceforge.openbabelJNI.OBMol_SetAromaticCorrected(swigCPtr, this);
  }

  public void SetSpinMultiplicityAssigned() {
    net.sourceforge.openbabelJNI.OBMol_SetSpinMultiplicityAssigned(swigCPtr, this);
  }

  public void SetFlags(int flags) {
    net.sourceforge.openbabelJNI.OBMol_SetFlags(swigCPtr, this, flags);
  }

  public void UnsetAromaticPerceived() {
    net.sourceforge.openbabelJNI.OBMol_UnsetAromaticPerceived(swigCPtr, this);
  }

  public void UnsetPartialChargesPerceived() {
    net.sourceforge.openbabelJNI.OBMol_UnsetPartialChargesPerceived(swigCPtr, this);
  }

  public void UnsetImplicitValencePerceived() {
    net.sourceforge.openbabelJNI.OBMol_UnsetImplicitValencePerceived(swigCPtr, this);
  }

  public void UnsetFlag(int flag) {
    net.sourceforge.openbabelJNI.OBMol_UnsetFlag(swigCPtr, this, flag);
  }

  public OBBase DoTransformations(SWIGTYPE_p_std__mapTstd__string_std__string_t pOptions) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_DoTransformations(swigCPtr, this, SWIGTYPE_p_std__mapTstd__string_std__string_t.getCPtr(pOptions));
    return (cPtr == 0) ? null : new OBBase(cPtr, false);
  }

  public static String ClassDescription() {
    return net.sourceforge.openbabelJNI.OBMol_ClassDescription();
  }

  public boolean Clear() {
    return net.sourceforge.openbabelJNI.OBMol_Clear(swigCPtr, this);
  }

  public void RenumberAtoms(SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t arg0) {
    net.sourceforge.openbabelJNI.OBMol_RenumberAtoms(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t.getCPtr(arg0));
  }

  public void ToInertialFrame(int conf, SWIGTYPE_p_double rmat) {
    net.sourceforge.openbabelJNI.OBMol_ToInertialFrame__SWIG_0(swigCPtr, this, conf, SWIGTYPE_p_double.getCPtr(rmat));
  }

  public void ToInertialFrame() {
    net.sourceforge.openbabelJNI.OBMol_ToInertialFrame__SWIG_1(swigCPtr, this);
  }

  public void Translate(vector3 v) {
    net.sourceforge.openbabelJNI.OBMol_Translate__SWIG_0(swigCPtr, this, vector3.getCPtr(v), v);
  }

  public void Translate(vector3 v, int conf) {
    net.sourceforge.openbabelJNI.OBMol_Translate__SWIG_1(swigCPtr, this, vector3.getCPtr(v), v, conf);
  }

  public void Rotate(SWIGTYPE_p_a_3__double u) {
    net.sourceforge.openbabelJNI.OBMol_Rotate__SWIG_0(swigCPtr, this, SWIGTYPE_p_a_3__double.getCPtr(u));
  }

  public void Rotate(SWIGTYPE_p_double m) {
    net.sourceforge.openbabelJNI.OBMol_Rotate__SWIG_1(swigCPtr, this, SWIGTYPE_p_double.getCPtr(m));
  }

  public void Rotate(SWIGTYPE_p_double m, int nconf) {
    net.sourceforge.openbabelJNI.OBMol_Rotate__SWIG_2(swigCPtr, this, SWIGTYPE_p_double.getCPtr(m), nconf);
  }

  public void Center() {
    net.sourceforge.openbabelJNI.OBMol_Center__SWIG_0(swigCPtr, this);
  }

  public boolean Kekulize() {
    return net.sourceforge.openbabelJNI.OBMol_Kekulize(swigCPtr, this);
  }

  public boolean PerceiveKekuleBonds() {
    return net.sourceforge.openbabelJNI.OBMol_PerceiveKekuleBonds(swigCPtr, this);
  }

  public void NewPerceiveKekuleBonds() {
    net.sourceforge.openbabelJNI.OBMol_NewPerceiveKekuleBonds(swigCPtr, this);
  }

  public boolean DeleteHydrogen(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_DeleteHydrogen(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean DeleteHydrogens() {
    return net.sourceforge.openbabelJNI.OBMol_DeleteHydrogens__SWIG_0(swigCPtr, this);
  }

  public boolean DeleteHydrogens(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_DeleteHydrogens__SWIG_1(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean DeleteNonPolarHydrogens() {
    return net.sourceforge.openbabelJNI.OBMol_DeleteNonPolarHydrogens(swigCPtr, this);
  }

  public boolean AddHydrogens(boolean polaronly, boolean correctForPH) {
    return net.sourceforge.openbabelJNI.OBMol_AddHydrogens__SWIG_0(swigCPtr, this, polaronly, correctForPH);
  }

  public boolean AddHydrogens(boolean polaronly) {
    return net.sourceforge.openbabelJNI.OBMol_AddHydrogens__SWIG_1(swigCPtr, this, polaronly);
  }

  public boolean AddHydrogens() {
    return net.sourceforge.openbabelJNI.OBMol_AddHydrogens__SWIG_2(swigCPtr, this);
  }

  public boolean AddHydrogens(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMol_AddHydrogens__SWIG_3(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean AddPolarHydrogens() {
    return net.sourceforge.openbabelJNI.OBMol_AddPolarHydrogens(swigCPtr, this);
  }

  public boolean StripSalts() {
    return net.sourceforge.openbabelJNI.OBMol_StripSalts(swigCPtr, this);
  }

  public boolean ConvertDativeBonds() {
    return net.sourceforge.openbabelJNI.OBMol_ConvertDativeBonds(swigCPtr, this);
  }

  public boolean CorrectForPH() {
    return net.sourceforge.openbabelJNI.OBMol_CorrectForPH(swigCPtr, this);
  }

  public boolean AssignSpinMultiplicity() {
    return net.sourceforge.openbabelJNI.OBMol_AssignSpinMultiplicity(swigCPtr, this);
  }

  public vector3 Center(int nconf) {
    return new vector3(net.sourceforge.openbabelJNI.OBMol_Center__SWIG_1(swigCPtr, this, nconf), true);
  }

  public void SetTorsion(OBAtom arg0, OBAtom arg1, OBAtom arg2, OBAtom arg3, double arg4) {
    net.sourceforge.openbabelJNI.OBMol_SetTorsion(swigCPtr, this, OBAtom.getCPtr(arg0), arg0, OBAtom.getCPtr(arg1), arg1, OBAtom.getCPtr(arg2), arg2, OBAtom.getCPtr(arg3), arg3, arg4);
  }

  public void FindSSSR() {
    net.sourceforge.openbabelJNI.OBMol_FindSSSR(swigCPtr, this);
  }

  public void FindRingAtomsAndBonds() {
    net.sourceforge.openbabelJNI.OBMol_FindRingAtomsAndBonds(swigCPtr, this);
  }

  public void FindChiralCenters() {
    net.sourceforge.openbabelJNI.OBMol_FindChiralCenters(swigCPtr, this);
  }

  public void FindChildren(vectorInt arg0, int arg1, int arg2) {
    net.sourceforge.openbabelJNI.OBMol_FindChildren__SWIG_0(swigCPtr, this, vectorInt.getCPtr(arg0), arg0, arg1, arg2);
  }

  public void FindChildren(SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t arg0, OBAtom arg1, OBAtom arg2) {
    net.sourceforge.openbabelJNI.OBMol_FindChildren__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBAtom_p_t.getCPtr(arg0), OBAtom.getCPtr(arg1), arg1, OBAtom.getCPtr(arg2), arg2);
  }

  public void FindLargestFragment(OBBitVec arg0) {
    net.sourceforge.openbabelJNI.OBMol_FindLargestFragment(swigCPtr, this, OBBitVec.getCPtr(arg0), arg0);
  }

  public void ContigFragList(vvInt arg0) {
    net.sourceforge.openbabelJNI.OBMol_ContigFragList(swigCPtr, this, vvInt.getCPtr(arg0), arg0);
  }

  public void Align(OBAtom arg0, OBAtom arg1, vector3 arg2, vector3 arg3) {
    net.sourceforge.openbabelJNI.OBMol_Align(swigCPtr, this, OBAtom.getCPtr(arg0), arg0, OBAtom.getCPtr(arg1), arg1, vector3.getCPtr(arg2), arg2, vector3.getCPtr(arg3), arg3);
  }

  public void ConnectTheDots() {
    net.sourceforge.openbabelJNI.OBMol_ConnectTheDots(swigCPtr, this);
  }

  public void PerceiveBondOrders() {
    net.sourceforge.openbabelJNI.OBMol_PerceiveBondOrders(swigCPtr, this);
  }

  public void FindTorsions() {
    net.sourceforge.openbabelJNI.OBMol_FindTorsions(swigCPtr, this);
  }

  public boolean GetGTDVector(vectorInt arg0) {
    return net.sourceforge.openbabelJNI.OBMol_GetGTDVector(swigCPtr, this, vectorInt.getCPtr(arg0), arg0);
  }

  public void GetGIVector(SWIGTYPE_p_std__vectorTunsigned_int_t arg0) {
    net.sourceforge.openbabelJNI.OBMol_GetGIVector(swigCPtr, this, SWIGTYPE_p_std__vectorTunsigned_int_t.getCPtr(arg0));
  }

  public void GetGIDVector(SWIGTYPE_p_std__vectorTunsigned_int_t arg0) {
    net.sourceforge.openbabelJNI.OBMol_GetGIDVector(swigCPtr, this, SWIGTYPE_p_std__vectorTunsigned_int_t.getCPtr(arg0));
  }

  public boolean Has2D() {
    return net.sourceforge.openbabelJNI.OBMol_Has2D(swigCPtr, this);
  }

  public boolean Has3D() {
    return net.sourceforge.openbabelJNI.OBMol_Has3D(swigCPtr, this);
  }

  public boolean HasNonZeroCoords() {
    return net.sourceforge.openbabelJNI.OBMol_HasNonZeroCoords(swigCPtr, this);
  }

  public boolean HasAromaticPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasAromaticPerceived(swigCPtr, this);
  }

  public boolean HasSSSRPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasSSSRPerceived(swigCPtr, this);
  }

  public boolean HasRingAtomsAndBondsPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasRingAtomsAndBondsPerceived(swigCPtr, this);
  }

  public boolean HasAtomTypesPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasAtomTypesPerceived(swigCPtr, this);
  }

  public boolean HasChiralityPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasChiralityPerceived(swigCPtr, this);
  }

  public boolean HasPartialChargesPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasPartialChargesPerceived(swigCPtr, this);
  }

  public boolean HasHybridizationPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasHybridizationPerceived(swigCPtr, this);
  }

  public boolean HasImplicitValencePerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasImplicitValencePerceived(swigCPtr, this);
  }

  public boolean HasKekulePerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasKekulePerceived(swigCPtr, this);
  }

  public boolean HasClosureBondsPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasClosureBondsPerceived(swigCPtr, this);
  }

  public boolean HasChainsPerceived() {
    return net.sourceforge.openbabelJNI.OBMol_HasChainsPerceived(swigCPtr, this);
  }

  public boolean HasHydrogensAdded() {
    return net.sourceforge.openbabelJNI.OBMol_HasHydrogensAdded(swigCPtr, this);
  }

  public boolean HasAromaticCorrected() {
    return net.sourceforge.openbabelJNI.OBMol_HasAromaticCorrected(swigCPtr, this);
  }

  public boolean IsCorrectedForPH() {
    return net.sourceforge.openbabelJNI.OBMol_IsCorrectedForPH(swigCPtr, this);
  }

  public boolean HasSpinMultiplicityAssigned() {
    return net.sourceforge.openbabelJNI.OBMol_HasSpinMultiplicityAssigned(swigCPtr, this);
  }

  public boolean IsChiral() {
    return net.sourceforge.openbabelJNI.OBMol_IsChiral(swigCPtr, this);
  }

  public boolean Empty() {
    return net.sourceforge.openbabelJNI.OBMol_Empty(swigCPtr, this);
  }

  public int NumConformers() {
    return net.sourceforge.openbabelJNI.OBMol_NumConformers(swigCPtr, this);
  }

  public void SetConformers(SWIGTYPE_p_std__vectorTdouble_p_t v) {
    net.sourceforge.openbabelJNI.OBMol_SetConformers(swigCPtr, this, SWIGTYPE_p_std__vectorTdouble_p_t.getCPtr(v));
  }

  public void AddConformer(SWIGTYPE_p_double f) {
    net.sourceforge.openbabelJNI.OBMol_AddConformer(swigCPtr, this, SWIGTYPE_p_double.getCPtr(f));
  }

  public void SetConformer(int i) {
    net.sourceforge.openbabelJNI.OBMol_SetConformer(swigCPtr, this, i);
  }

  public void CopyConformer(SWIGTYPE_p_double arg0, int arg1) {
    net.sourceforge.openbabelJNI.OBMol_CopyConformer(swigCPtr, this, SWIGTYPE_p_double.getCPtr(arg0), arg1);
  }

  public void DeleteConformer(int arg0) {
    net.sourceforge.openbabelJNI.OBMol_DeleteConformer(swigCPtr, this, arg0);
  }

  public SWIGTYPE_p_double GetConformer(int i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_GetConformer(swigCPtr, this, i);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public SWIGTYPE_p_double BeginConformer(SWIGTYPE_p_std__vectorTdouble_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_BeginConformer(swigCPtr, this, SWIGTYPE_p_std__vectorTdouble_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public SWIGTYPE_p_double NextConformer(SWIGTYPE_p_std__vectorTdouble_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NextConformer(swigCPtr, this, SWIGTYPE_p_std__vectorTdouble_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTdouble_p_t GetConformers() {
    return new SWIGTYPE_p_std__vectorTdouble_p_t(net.sourceforge.openbabelJNI.OBMol_GetConformers(swigCPtr, this), false);
  }

  public OBAtom BeginAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_BeginAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom NextAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NextAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBNodeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBBond BeginBond(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_BeginBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond NextBond(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NextBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBResidue BeginResidue(SWIGTYPE_p_std__vectorTOpenBabel__OBResidue_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_BeginResidue(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBResidue_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public OBResidue NextResidue(SWIGTYPE_p_std__vectorTOpenBabel__OBResidue_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NextResidue(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBResidue_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public SWIGTYPE_p_OpenBabel__OBInternalCoord BeginInternalCoord(SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_BeginInternalCoord(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new SWIGTYPE_p_OpenBabel__OBInternalCoord(cPtr, false);
  }

  public SWIGTYPE_p_OpenBabel__OBInternalCoord NextInternalCoord(SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMol_NextInternalCoord(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBInternalCoord_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new SWIGTYPE_p_OpenBabel__OBInternalCoord(cPtr, false);
  }

}
