/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class OBMolTorsionIter {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected OBMolTorsionIter(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(OBMolTorsionIter obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_OBMolTorsionIter(swigCPtr);
    }
    swigCPtr = 0;
  }

  public OBMolTorsionIter() {
    this(openbabelJNI.new_OBMolTorsionIter__SWIG_0(), true);
  }

  public OBMolTorsionIter(OBMol mol) {
    this(openbabelJNI.new_OBMolTorsionIter__SWIG_1(OBMol.getCPtr(mol), mol), true);
  }

  public OBMolTorsionIter(OBMolTorsionIter ai) {
    this(openbabelJNI.new_OBMolTorsionIter__SWIG_2(OBMolTorsionIter.getCPtr(ai), ai), true);
  }

  public boolean good() {
    return openbabelJNI.OBMolTorsionIter_good(swigCPtr, this);
  }

  public OBMolTorsionIter inc() {
    return new OBMolTorsionIter(openbabelJNI.OBMolTorsionIter_inc(swigCPtr, this), false);
  }

  public vectorUnsignedInt __ref__() {
    return new vectorUnsignedInt(openbabelJNI.OBMolTorsionIter___ref__(swigCPtr, this), true);
  }

}
