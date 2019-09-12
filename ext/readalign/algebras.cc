//
// algebras.cc - extended real types
//
// Gerton Lunter, 27/8/04
//
//


#include "algebras.h"


BFMantissa *BFloat::aConversionLookup;           // Actual location of the static members of BFloat class
double *BFloat::aDoubleConversionLookup;


_BFloatInitialize _dummyInitializer;             // This initializes aConversionLookup and aDoubleConversionLookup


_BFloatInitialize::_BFloatInitialize() {

  BFloat::aConversionLookup = new BFMantissa[cBFloatConvTableSize];
  BFloat::aDoubleConversionLookup = new double[cBFloatDoubleConvTableSize];

  BFMantissa iBFM = 1.0;
  for (int i = 0; i < cBFloatConvTableSize; i++) {
    BFloat::aConversionLookup[ i ] = iBFM;
    iBFM *= cBFloatRangeInv;
  }

  for (int i = 0; i < cBFloatDoubleConvTableSize; i++) {
    BFloat::aDoubleConversionLookup[ i ] = exp( (i-cBFloatDoubleConvTableSize/2) * logcBFloatRange );
  }

}
