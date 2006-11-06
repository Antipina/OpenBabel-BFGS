/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolAtomIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolAtomIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolAtomIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_OBMolAtomIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolAtomIter() {
    this(net.sourceforge.openbabelJNI.new_OBMolAtomIter__SWIG_0(), true);
  }

  public OBMolAtomIter(OBMol mol) {
    this(net.sourceforge.openbabelJNI.new_OBMolAtomIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolAtomIter(OBMolAtomIter ai) {
    this(net.sourceforge.openbabelJNI.new_OBMolAtomIter__SWIG_3(OBMolAtomIter.getCPtr(ai), ai), true);
  }

  public boolean good() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_good(swigCPtr, this);
  }

  public OBMolAtomIter inc(int arg0) {
    return new OBMolAtomIter(net.sourceforge.openbabelJNI.OBMolAtomIter_inc(swigCPtr, this, arg0), true);
  }

  public OBAtom deref() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_deref(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom __ref__() {
    return new OBAtom(net.sourceforge.openbabelJNI.OBMolAtomIter___ref__(swigCPtr, this), false);
  }

  public void Clear() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_Clear(swigCPtr, this);
  }

  public void SetIdx(int idx) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetIdx(swigCPtr, this, idx);
  }

  public void SetHyb(int hyb) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetHyb(swigCPtr, this, hyb);
  }

  public void SetAtomicNum(int atomicnum) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetAtomicNum(swigCPtr, this, atomicnum);
  }

  public void SetIsotope(long iso) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetIsotope(swigCPtr, this, iso);
  }

  public void SetImplicitValence(int val) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetImplicitValence(swigCPtr, this, val);
  }

  public void IncrementImplicitValence() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_IncrementImplicitValence(swigCPtr, this);
  }

  public void DecrementImplicitValence() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_DecrementImplicitValence(swigCPtr, this);
  }

  public void SetFormalCharge(int fcharge) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetFormalCharge(swigCPtr, this, fcharge);
  }

  public void SetSpinMultiplicity(short spin) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetSpinMultiplicity(swigCPtr, this, spin);
  }

  public void SetType(String type) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetType__SWIG_0(swigCPtr, this, type);
  }

  public void SetType(SWIGTYPE_p_std__string type) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetType__SWIG_1(swigCPtr, this, SWIGTYPE_p_std__string.getCPtr(type));
  }

  public void SetPartialCharge(double pcharge) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetPartialCharge(swigCPtr, this, pcharge);
  }

  public void SetVector(vector3 v) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetVector__SWIG_0(swigCPtr, this, vector3.getCPtr(v), v);
  }

  public void SetVector(double x, double y, double z) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetVector__SWIG_1(swigCPtr, this, x, y, z);
  }

  public void SetVector() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetVector__SWIG_2(swigCPtr, this);
  }

  public void SetCoordPtr(SWIGTYPE_p_p_double c) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetCoordPtr(swigCPtr, this, SWIGTYPE_p_p_double.getCPtr(c));
  }

  public void SetResidue(OBResidue res) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetResidue(swigCPtr, this, OBResidue.getCPtr(res), res);
  }

  public void SetAromatic() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetAromatic(swigCPtr, this);
  }

  public void UnsetAromatic() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_UnsetAromatic(swigCPtr, this);
  }

  public void SetClockwiseStereo() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetClockwiseStereo(swigCPtr, this);
  }

  public void SetAntiClockwiseStereo() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetAntiClockwiseStereo(swigCPtr, this);
  }

  public void SetPositiveStereo() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetPositiveStereo(swigCPtr, this);
  }

  public void SetNegativeStereo() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetNegativeStereo(swigCPtr, this);
  }

  public void UnsetStereo() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_UnsetStereo(swigCPtr, this);
  }

  public void SetInRing() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetInRing(swigCPtr, this);
  }

  public void SetChiral() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetChiral(swigCPtr, this);
  }

  public void ClearCoordPtr() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_ClearCoordPtr(swigCPtr, this);
  }

  public int GetFormalCharge() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetFormalCharge(swigCPtr, this);
  }

  public long GetAtomicNum() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetAtomicNum(swigCPtr, this);
  }

  public int GetIsotope() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetIsotope(swigCPtr, this);
  }

  public int GetSpinMultiplicity() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetSpinMultiplicity(swigCPtr, this);
  }

  public double GetAtomicMass() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetAtomicMass(swigCPtr, this);
  }

  public double GetExactMass() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetExactMass(swigCPtr, this);
  }

  public long GetIdx() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetIdx(swigCPtr, this);
  }

  public long GetCoordinateIdx() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetCoordinateIdx(swigCPtr, this);
  }

  public long GetCIdx() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetCIdx(swigCPtr, this);
  }

  public long GetValence() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetValence(swigCPtr, this);
  }

  public long GetHyb() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetHyb(swigCPtr, this);
  }

  public long GetImplicitValence() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetImplicitValence(swigCPtr, this);
  }

  public long GetHvyValence() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetHvyValence(swigCPtr, this);
  }

  public long GetHeteroValence() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetHeteroValence(swigCPtr, this);
  }

  public String GetType() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetType(swigCPtr, this);
  }

  public double GetX() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetX(swigCPtr, this);
  }

  public double x() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_x(swigCPtr, this);
  }

  public double GetY() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetY(swigCPtr, this);
  }

  public double y() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_y(swigCPtr, this);
  }

  public double GetZ() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetZ(swigCPtr, this);
  }

  public double z() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_z(swigCPtr, this);
  }

  public SWIGTYPE_p_double GetCoordinate() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetCoordinate(swigCPtr, this);
    return (cPtr == 0) ? null : new SWIGTYPE_p_double(cPtr, false);
  }

  public vector3 GetVector() {
    return new vector3(net.sourceforge.openbabelJNI.OBMolAtomIter_GetVector(swigCPtr, this), false);
  }

  public double GetPartialCharge() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetPartialCharge(swigCPtr, this);
  }

  public OBResidue GetResidue() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetResidue(swigCPtr, this);
    return (cPtr == 0) ? null : new OBResidue(cPtr, false);
  }

  public boolean GetNewBondVector(vector3 v, double length) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetNewBondVector(swigCPtr, this, vector3.getCPtr(v), v, length);
  }

  public OBBond GetBond(OBAtom arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetBond(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom GetNextAtom() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetNextAtom(swigCPtr, this);
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator BeginBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator(net.sourceforge.openbabelJNI.OBMolAtomIter_BeginBonds(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator EndBonds() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator(net.sourceforge.openbabelJNI.OBMolAtomIter_EndBonds(swigCPtr, this), true);
  }

  public OBBond BeginBond(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_BeginBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBBond NextBond(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator i) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_NextBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(i));
    return (cPtr == 0) ? null : new OBBond(cPtr, false);
  }

  public OBAtom BeginNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_BeginNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public OBAtom NextNbrAtom(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_NextNbrAtom(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBAtom(cPtr, false);
  }

  public double GetDistance(int index) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetDistance__SWIG_0(swigCPtr, this, index);
  }

  public double GetDistance(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetDistance__SWIG_1(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public double GetAngle(int b, int c) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetAngle__SWIG_0(swigCPtr, this, b, c);
  }

  public double GetAngle(OBAtom b, OBAtom c) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_GetAngle__SWIG_1(swigCPtr, this, OBAtom.getCPtr(b), b, OBAtom.getCPtr(c), c);
  }

  public void NewResidue() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_NewResidue(swigCPtr, this);
  }

  public void DeleteResidue() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_DeleteResidue(swigCPtr, this);
  }

  public void AddBond(OBBond bond) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_AddBond(swigCPtr, this, OBBond.getCPtr(bond), bond);
  }

  public void InsertBond(SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator i, OBBond bond) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_InsertBond(swigCPtr, this, SWIGTYPE_p_std__vectorTOpenBabel__OBEdgeBase_p_t__iterator.getCPtr(i), OBBond.getCPtr(bond), bond);
  }

  public boolean DeleteBond(OBBond arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_DeleteBond(swigCPtr, this, OBBond.getCPtr(arg0), arg0);
  }

  public void ClearBond() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_ClearBond(swigCPtr, this);
  }

  public long CountFreeOxygens() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_CountFreeOxygens(swigCPtr, this);
  }

  public long ImplicitHydrogenCount() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_ImplicitHydrogenCount(swigCPtr, this);
  }

  public long ExplicitHydrogenCount(boolean ExcludeIsotopes) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_ExplicitHydrogenCount__SWIG_0(swigCPtr, this, ExcludeIsotopes);
  }

  public long ExplicitHydrogenCount() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_ExplicitHydrogenCount__SWIG_1(swigCPtr, this);
  }

  public long MemberOfRingCount() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_MemberOfRingCount(swigCPtr, this);
  }

  public long MemberOfRingSize() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_MemberOfRingSize(swigCPtr, this);
  }

  public long CountRingBonds() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_CountRingBonds(swigCPtr, this);
  }

  public double SmallestBondAngle() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_SmallestBondAngle(swigCPtr, this);
  }

  public double AverageBondAngle() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_AverageBondAngle(swigCPtr, this);
  }

  public long BOSum() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_BOSum(swigCPtr, this);
  }

  public long KBOSum() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_KBOSum(swigCPtr, this);
  }

  public boolean HtoMethyl() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HtoMethyl(swigCPtr, this);
  }

  public boolean SetHybAndGeom(int arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_SetHybAndGeom(swigCPtr, this, arg0);
  }

  public void ForceNoH() {
    net.sourceforge.openbabelJNI.OBMolAtomIter_ForceNoH(swigCPtr, this);
  }

  public boolean HasNoHForced() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasNoHForced(swigCPtr, this);
  }

  public boolean HasResidue() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasResidue(swigCPtr, this);
  }

  public boolean IsHydrogen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsHydrogen(swigCPtr, this);
  }

  public boolean IsCarbon() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsCarbon(swigCPtr, this);
  }

  public boolean IsNitrogen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsNitrogen(swigCPtr, this);
  }

  public boolean IsOxygen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsOxygen(swigCPtr, this);
  }

  public boolean IsSulfur() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsSulfur(swigCPtr, this);
  }

  public boolean IsPhosphorus() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsPhosphorus(swigCPtr, this);
  }

  public boolean IsAromatic() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsAromatic(swigCPtr, this);
  }

  public boolean IsInRing() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsInRing(swigCPtr, this);
  }

  public boolean IsInRingSize(int arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsInRingSize(swigCPtr, this, arg0);
  }

  public boolean IsHeteroatom() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsHeteroatom(swigCPtr, this);
  }

  public boolean IsNotCorH() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsNotCorH(swigCPtr, this);
  }

  public boolean IsConnected(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsConnected(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneThree(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsOneThree(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsOneFour(OBAtom arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsOneFour(swigCPtr, this, OBAtom.getCPtr(arg0), arg0);
  }

  public boolean IsCarboxylOxygen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsCarboxylOxygen(swigCPtr, this);
  }

  public boolean IsPhosphateOxygen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsPhosphateOxygen(swigCPtr, this);
  }

  public boolean IsSulfateOxygen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsSulfateOxygen(swigCPtr, this);
  }

  public boolean IsNitroOxygen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsNitroOxygen(swigCPtr, this);
  }

  public boolean IsAmideNitrogen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsAmideNitrogen(swigCPtr, this);
  }

  public boolean IsPolarHydrogen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsPolarHydrogen(swigCPtr, this);
  }

  public boolean IsNonPolarHydrogen() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsNonPolarHydrogen(swigCPtr, this);
  }

  public boolean IsAromaticNOxide() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsAromaticNOxide(swigCPtr, this);
  }

  public boolean IsChiral() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsChiral(swigCPtr, this);
  }

  public boolean IsAxial() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsAxial(swigCPtr, this);
  }

  public boolean IsClockwise() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsClockwise(swigCPtr, this);
  }

  public boolean IsAntiClockwise() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsAntiClockwise(swigCPtr, this);
  }

  public boolean IsPositiveStereo() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsPositiveStereo(swigCPtr, this);
  }

  public boolean IsNegativeStereo() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsNegativeStereo(swigCPtr, this);
  }

  public boolean HasChiralitySpecified() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasChiralitySpecified(swigCPtr, this);
  }

  public boolean HasChiralVolume() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasChiralVolume(swigCPtr, this);
  }

  public boolean IsHbondAcceptor() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsHbondAcceptor(swigCPtr, this);
  }

  public boolean IsHbondDonor() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsHbondDonor(swigCPtr, this);
  }

  public boolean IsHbondDonorH() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_IsHbondDonorH(swigCPtr, this);
  }

  public boolean HasAlphaBetaUnsat(boolean includePandS) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasAlphaBetaUnsat__SWIG_0(swigCPtr, this, includePandS);
  }

  public boolean HasAlphaBetaUnsat() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasAlphaBetaUnsat__SWIG_1(swigCPtr, this);
  }

  public boolean HasBondOfOrder(long arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasBondOfOrder(swigCPtr, this, arg0);
  }

  public int CountBondsOfOrder(long arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_CountBondsOfOrder(swigCPtr, this, arg0);
  }

  public boolean HasNonSingleBond() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasNonSingleBond(swigCPtr, this);
  }

  public boolean HasSingleBond() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasSingleBond(swigCPtr, this);
  }

  public boolean HasDoubleBond() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasDoubleBond(swigCPtr, this);
  }

  public boolean HasAromaticBond() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasAromaticBond(swigCPtr, this);
  }

  public boolean MatchesSMARTS(String arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_MatchesSMARTS(swigCPtr, this, arg0);
  }

  public void setVisit(boolean value) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_Visit_set(swigCPtr, this, value);
  }

  public boolean getVisit() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_Visit_get(swigCPtr, this);
  }

  public OBGraphBase GetParent() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetParent(swigCPtr, this);
    return (cPtr == 0) ? null : new OBGraphBase(cPtr, false);
  }

  public void SetParent(OBGraphBase arg0) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetParent(swigCPtr, this, OBGraphBase.getCPtr(arg0), arg0);
  }

  public void AddEdge(OBEdgeBase b) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_AddEdge(swigCPtr, this, OBEdgeBase.getCPtr(b), b);
  }

  public void Error(int f) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_Error(swigCPtr, this, f);
  }

  public void SetMatch(OBNodeBase arg0) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetMatch(swigCPtr, this, OBNodeBase.getCPtr(arg0), arg0);
  }

  public boolean Eval(OBNodeBase arg0) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_Eval(swigCPtr, this, OBNodeBase.getCPtr(arg0), arg0);
  }

  public OBNodeBase GetMatch() {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetMatch(swigCPtr, this);
    return (cPtr == 0) ? null : new OBNodeBase(cPtr, false);
  }

  public OBBase DoTransformations(SWIGTYPE_p_std__mapTstd__string_std__string_t arg0) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_DoTransformations(swigCPtr, this, SWIGTYPE_p_std__mapTstd__string_std__string_t.getCPtr(arg0));
    return (cPtr == 0) ? null : new OBBase(cPtr, false);
  }

  public String ClassDescription() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_ClassDescription(swigCPtr, this);
  }

  public boolean HasData(long type) {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_HasData__SWIG_2(swigCPtr, this, type);
  }

  public void DeleteData(long type) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_DeleteData__SWIG_0(swigCPtr, this, type);
  }

  public void DeleteData(OBGenericData arg0) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_DeleteData__SWIG_1(swigCPtr, this, OBGenericData.getCPtr(arg0), arg0);
  }

  public void DeleteData(vectorData arg0) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_DeleteData__SWIG_2(swigCPtr, this, vectorData.getCPtr(arg0), arg0);
  }

  public void SetData(OBGenericData d) {
    net.sourceforge.openbabelJNI.OBMolAtomIter_SetData(swigCPtr, this, OBGenericData.getCPtr(d), d);
  }

  public long DataSize() {
    return net.sourceforge.openbabelJNI.OBMolAtomIter_DataSize(swigCPtr, this);
  }

  public OBGenericData GetData(long type) {
    long cPtr = net.sourceforge.openbabelJNI.OBMolAtomIter_GetData__SWIG_0(swigCPtr, this, type);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public vectorData GetData() {
    return new vectorData(net.sourceforge.openbabelJNI.OBMolAtomIter_GetData__SWIG_3(swigCPtr, this), false);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator BeginData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(net.sourceforge.openbabelJNI.OBMolAtomIter_BeginData(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator EndData() {
    return new SWIGTYPE_p_std__vectorTOpenBabel__OBGenericData_p_t__iterator(net.sourceforge.openbabelJNI.OBMolAtomIter_EndData(swigCPtr, this), true);
  }

}
