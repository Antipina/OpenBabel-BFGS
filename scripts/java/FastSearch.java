/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.30
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class FastSearch {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected FastSearch(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(FastSearch obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_FastSearch(swigCPtr);
    }
    swigCPtr = 0;
  }

  public String ReadIndex(SWIGTYPE_p_std__istream pIndexstream) {
    return openbabelJNI.FastSearch_ReadIndex(swigCPtr, this, SWIGTYPE_p_std__istream.getCPtr(pIndexstream));
  }

  public boolean Find(OBBase pOb, SWIGTYPE_p_std__vectorTunsigned_int_t SeekPositions, long MaxCandidates) {
    return openbabelJNI.FastSearch_Find(swigCPtr, this, OBBase.getCPtr(pOb), pOb, SWIGTYPE_p_std__vectorTunsigned_int_t.getCPtr(SeekPositions), MaxCandidates);
  }

  public boolean FindSimilar(OBBase pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t SeekposMap, double MinTani) {
    return openbabelJNI.FastSearch_FindSimilar__SWIG_0(swigCPtr, this, OBBase.getCPtr(pOb), pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t.getCPtr(SeekposMap), MinTani);
  }

  public boolean FindSimilar(OBBase pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t SeekposMap, int nCandidates) {
    return openbabelJNI.FastSearch_FindSimilar__SWIG_1(swigCPtr, this, OBBase.getCPtr(pOb), pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t.getCPtr(SeekposMap), nCandidates);
  }

  public boolean FindSimilar(OBBase pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t SeekposMap) {
    return openbabelJNI.FastSearch_FindSimilar__SWIG_2(swigCPtr, this, OBBase.getCPtr(pOb), pOb, SWIGTYPE_p_std__multimapTdouble_unsigned_int_t.getCPtr(SeekposMap));
  }

  public OBFingerprint GetFingerprint() {
    long cPtr = openbabelJNI.FastSearch_GetFingerprint(swigCPtr, this);
    return (cPtr == 0) ? null : new OBFingerprint(cPtr, false);
  }

  public FptIndexHeader GetIndexHeader() {
    return new FptIndexHeader(openbabelJNI.FastSearch_GetIndexHeader(swigCPtr, this), false);
  }

  public FastSearch() {
    this(openbabelJNI.new_FastSearch(), true);
  }

}
