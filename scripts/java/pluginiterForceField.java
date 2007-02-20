/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.31
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class pluginiterForceField {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected pluginiterForceField(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(pluginiterForceField obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      openbabelJNI.delete_pluginiterForceField(swigCPtr);
    }
    swigCPtr = 0;
  }

  public void Register(OBForceField pType, String ID, boolean IsDefault) {
    openbabelJNI.pluginiterForceField_Register(swigCPtr, this, OBForceField.getCPtr(pType), pType, ID, IsDefault);
  }

  public OBForceField FindType(String ID) {
    long cPtr = openbabelJNI.pluginiterForceField_FindType(swigCPtr, this, ID);
    return (cPtr == 0) ? null : new OBForceField(cPtr, false);
  }

  public OBForceField FindDefaultType() {
    long cPtr = openbabelJNI.pluginiterForceField_FindDefaultType(swigCPtr, this);
    return (cPtr == 0) ? null : new OBForceField(cPtr, false);
  }

  public String ID() {
    return openbabelJNI.pluginiterForceField_ID(swigCPtr, this);
  }

  public String Description() {
    return openbabelJNI.pluginiterForceField_Description(swigCPtr, this);
  }

  public void ToStart() {
    openbabelJNI.pluginiterForceField_ToStart(swigCPtr, this);
  }

  public OBForceField __deref__() {
    long cPtr = openbabelJNI.pluginiterForceField___deref__(swigCPtr, this);
    return (cPtr == 0) ? null : new OBForceField(cPtr, false);
  }

  public OBForceField __ref__() {
    return new OBForceField(openbabelJNI.pluginiterForceField___ref__(swigCPtr, this), false);
  }

  public pluginiterForceField() {
    this(openbabelJNI.new_pluginiterForceField(), true);
  }

  public pluginiterForceField Iter() {
    return new pluginiterForceField(openbabelJNI.pluginiterForceField_Iter(swigCPtr, this), false);
  }

  public OBForceField FindForceField(String ID) {
    long cPtr = openbabelJNI.pluginiterForceField_FindForceField__SWIG_0(swigCPtr, this, ID);
    return (cPtr == 0) ? null : new OBForceField(cPtr, false);
  }

  public String GetUnit() {
    return openbabelJNI.pluginiterForceField_GetUnit(swigCPtr, this);
  }

  public boolean Setup(OBMol mol) {
    return openbabelJNI.pluginiterForceField_Setup(swigCPtr, this, OBMol.getCPtr(mol), mol);
  }

  public void UpdateCoordinates(OBMol mol) {
    openbabelJNI.pluginiterForceField_UpdateCoordinates(swigCPtr, this, OBMol.getCPtr(mol), mol);
  }

  public double Energy() {
    return openbabelJNI.pluginiterForceField_Energy(swigCPtr, this);
  }

  public double E_Bond() {
    return openbabelJNI.pluginiterForceField_E_Bond(swigCPtr, this);
  }

  public double E_Angle() {
    return openbabelJNI.pluginiterForceField_E_Angle(swigCPtr, this);
  }

  public double E_StrBnd() {
    return openbabelJNI.pluginiterForceField_E_StrBnd(swigCPtr, this);
  }

  public double E_Torsion() {
    return openbabelJNI.pluginiterForceField_E_Torsion(swigCPtr, this);
  }

  public double E_OOP() {
    return openbabelJNI.pluginiterForceField_E_OOP(swigCPtr, this);
  }

  public double E_VDW() {
    return openbabelJNI.pluginiterForceField_E_VDW(swigCPtr, this);
  }

  public double E_Electrostatic() {
    return openbabelJNI.pluginiterForceField_E_Electrostatic(swigCPtr, this);
  }

  public boolean SetLogFile(SWIGTYPE_p_std__ostream pos) {
    return openbabelJNI.pluginiterForceField_SetLogFile(swigCPtr, this, SWIGTYPE_p_std__ostream.getCPtr(pos));
  }

  public boolean SetLogLevel(int level) {
    return openbabelJNI.pluginiterForceField_SetLogLevel(swigCPtr, this, level);
  }

  public int GetLogLevel() {
    return openbabelJNI.pluginiterForceField_GetLogLevel(swigCPtr, this);
  }

  public void DistanceGeometry() {
    openbabelJNI.pluginiterForceField_DistanceGeometry(swigCPtr, this);
  }

  public void GenerateCoordinates() {
    openbabelJNI.pluginiterForceField_GenerateCoordinates(swigCPtr, this);
  }

  public void SystematicRotorSearch() {
    openbabelJNI.pluginiterForceField_SystematicRotorSearch(swigCPtr, this);
  }

  public vector3 LineSearch(OBAtom atom, vector3 direction) {
    return new vector3(openbabelJNI.pluginiterForceField_LineSearch(swigCPtr, this, OBAtom.getCPtr(atom), atom, vector3.getCPtr(direction), direction), true);
  }

  public void SteepestDescent(int steps, int method) {
    openbabelJNI.pluginiterForceField_SteepestDescent__SWIG_0(swigCPtr, this, steps, method);
  }

  public void SteepestDescent(int steps) {
    openbabelJNI.pluginiterForceField_SteepestDescent__SWIG_1(swigCPtr, this, steps);
  }

  public void ConjugateGradients(int steps, int method) {
    openbabelJNI.pluginiterForceField_ConjugateGradients__SWIG_0(swigCPtr, this, steps, method);
  }

  public void ConjugateGradients(int steps) {
    openbabelJNI.pluginiterForceField_ConjugateGradients__SWIG_1(swigCPtr, this, steps);
  }

  public vector3 ValidateLineSearch(OBAtom atom, vector3 direction) {
    return new vector3(openbabelJNI.pluginiterForceField_ValidateLineSearch(swigCPtr, this, OBAtom.getCPtr(atom), atom, vector3.getCPtr(direction), direction), true);
  }

  public void ValidateSteepestDescent(int steps) {
    openbabelJNI.pluginiterForceField_ValidateSteepestDescent(swigCPtr, this, steps);
  }

  public void ValidateConjugateGradients(int steps) {
    openbabelJNI.pluginiterForceField_ValidateConjugateGradients(swigCPtr, this, steps);
  }

  public boolean Validate() {
    return openbabelJNI.pluginiterForceField_Validate(swigCPtr, this);
  }

  public boolean ValidateGradients() {
    return openbabelJNI.pluginiterForceField_ValidateGradients(swigCPtr, this);
  }

  public vector3 ValidateGradientError(vector3 numgrad, vector3 anagrad) {
    return new vector3(openbabelJNI.pluginiterForceField_ValidateGradientError(swigCPtr, this, vector3.getCPtr(numgrad), numgrad, vector3.getCPtr(anagrad), anagrad), true);
  }

  public double VectorLengthDerivative(vector3 a, vector3 b) {
    return openbabelJNI.pluginiterForceField_VectorLengthDerivative(swigCPtr, this, vector3.getCPtr(a), a, vector3.getCPtr(b), b);
  }

  public double VectorAngleDerivative(vector3 a, vector3 b, vector3 c) {
    return openbabelJNI.pluginiterForceField_VectorAngleDerivative(swigCPtr, this, vector3.getCPtr(a), a, vector3.getCPtr(b), b, vector3.getCPtr(c), c);
  }

  public double VectorTorsionDerivative(vector3 a, vector3 b, vector3 c, vector3 d) {
    return openbabelJNI.pluginiterForceField_VectorTorsionDerivative(swigCPtr, this, vector3.getCPtr(a), a, vector3.getCPtr(b), b, vector3.getCPtr(c), c, vector3.getCPtr(d), d);
  }

}