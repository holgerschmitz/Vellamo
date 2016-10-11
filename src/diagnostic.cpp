/*
 * diagnostic.cpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"

FieldDiagnostic::IndexType FieldDiagnostic::getGlobalMin()
{
  return IndexType(0);
}

FieldDiagnostic::IndexType FieldDiagnostic::getGlobalMax()
{
  return Vellamo::getGlobalMax();
}


