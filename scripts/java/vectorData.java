/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.29
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class vectorData {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected vectorData(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(vectorData obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      net.sourceforge.openbabelJNI.delete_vectorData(swigCPtr);
    }
    swigCPtr = 0;
  }

  public vectorData() {
    this(net.sourceforge.openbabelJNI.new_vectorData__SWIG_0(), true);
  }

  public vectorData(long n) {
    this(net.sourceforge.openbabelJNI.new_vectorData__SWIG_1(n), true);
  }

  public long size() {
    return net.sourceforge.openbabelJNI.vectorData_size(swigCPtr);
  }

  public long capacity() {
    return net.sourceforge.openbabelJNI.vectorData_capacity(swigCPtr);
  }

  public void reserve(long n) {
    net.sourceforge.openbabelJNI.vectorData_reserve(swigCPtr, n);
  }

  public boolean isEmpty() {
    return net.sourceforge.openbabelJNI.vectorData_isEmpty(swigCPtr);
  }

  public void clear() {
    net.sourceforge.openbabelJNI.vectorData_clear(swigCPtr);
  }

  public void add(OBGenericData x) {
    net.sourceforge.openbabelJNI.vectorData_add(swigCPtr, OBGenericData.getCPtr(x));
  }

  public OBGenericData get(int i) {
    long cPtr = net.sourceforge.openbabelJNI.vectorData_get(swigCPtr, i);
    return (cPtr == 0) ? null : new OBGenericData(cPtr, false);
  }

  public void set(int i, OBGenericData x) {
    net.sourceforge.openbabelJNI.vectorData_set(swigCPtr, i, OBGenericData.getCPtr(x));
  }

}
