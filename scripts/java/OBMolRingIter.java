/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolRingIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolRingIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolRingIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMolRingIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolRingIter() {
    this(openbabelJNI.new_OBMolRingIter__SWIG_0(), true);
  }

  public OBMolRingIter(OBMol mol) {
    this(openbabelJNI.new_OBMolRingIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolRingIter(OBMolRingIter ri) {
    this(openbabelJNI.new_OBMolRingIter__SWIG_3(OBMolRingIter.getCPtr(ri), ri), true);
  }

  public OBRing __deref__() {
    long cPtr = openbabelJNI.OBMolRingIter___deref__(swigCPtr, this);
    return (cPtr == 0) ? null : new OBRing(cPtr, false);
  }

  public OBRing __ref__() {
    return new OBRing(openbabelJNI.OBMolRingIter___ref__(swigCPtr, this), false);
  }

  public void set_path(vectorInt value) {
    openbabelJNI.OBMolRingIter__path_set(swigCPtr, this, vectorInt.getCPtr(value), value);
  }

  public vectorInt get_path() {
    return new vectorInt(openbabelJNI.OBMolRingIter__path_get(swigCPtr, this), false);
  }

  public void set_pathset(OBBitVec value) {
    openbabelJNI.OBMolRingIter__pathset_set(swigCPtr, this, OBBitVec.getCPtr(value), value);
  }

  public OBBitVec get_pathset() {
    return new OBBitVec(openbabelJNI.OBMolRingIter__pathset_get(swigCPtr, this), false);
  }

  public boolean findCenterAndNormal(vector3 center, vector3 norm1, vector3 norm2) {
    return openbabelJNI.OBMolRingIter_findCenterAndNormal(swigCPtr, this, vector3.getCPtr(center), center, vector3.getCPtr(norm1), norm1, vector3.getCPtr(norm2), norm2);
  }

  public int Size() {
    return openbabelJNI.OBMolRingIter_Size(swigCPtr, this);
  }

  public int PathSize() {
    return openbabelJNI.OBMolRingIter_PathSize(swigCPtr, this);
  }

  public boolean IsMember(OBAtom a) {
    return openbabelJNI.OBMolRingIter_IsMember__SWIG_0(swigCPtr, this, OBAtom.getCPtr(a), a);
  }

  public boolean IsMember(OBBond b) {
    return openbabelJNI.OBMolRingIter_IsMember__SWIG_1(swigCPtr, this, OBBond.getCPtr(b), b);
  }

  public boolean IsAromatic() {
    return openbabelJNI.OBMolRingIter_IsAromatic(swigCPtr, this);
  }

  public boolean IsInRing(int i) {
    return openbabelJNI.OBMolRingIter_IsInRing(swigCPtr, this, i);
  }

  public void SetParent(OBMol m) {
    openbabelJNI.OBMolRingIter_SetParent(swigCPtr, this, OBMol.getCPtr(m), m);
  }

  public OBMol GetParent() {
    long cPtr = openbabelJNI.OBMolRingIter_GetParent(swigCPtr, this);
    return (cPtr == 0) ? null : new OBMol(cPtr, false);
  }

}
