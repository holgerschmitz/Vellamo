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

schnek::pHdfAttributes FieldDiagnostic::getAttributes()
{
  static double time;
  static Vector dx;
  static Index grid_size;

  time = Vellamo::getT();
  dx = Vellamo::getDx();
  grid_size = Vellamo::getGlobalMax();

  schnek::pHdfAttributes attributes(new schnek::HdfAttributes());
  attributes->set("time", time);
  attributes->set("dx", dx[0], DIMENSION);
  attributes->set("grid_size", grid_size[0], DIMENSION);
  return attributes;
}

