/*
 * boundary.cpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#include "boundary.hpp"

void BoundaryCondition::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addArrayParameter("low_",applyLo,Index(0));
  blockPars.addArrayParameter("high_",applyHi,Index(0));
}

void BoundaryCondition::apply(Field& f)
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  if (bool(applyLo[0]) && subdivision.isBoundLo(0)) applyLoX(f);
  if (bool(applyLo[1]) && subdivision.isBoundLo(1)) applyLoY(f);
  if (bool(applyHi[0]) && subdivision.isBoundHi(0)) applyHiX(f);
  if (bool(applyHi[1]) && subdivision.isBoundHi(1)) applyHiY(f);
}

void ZeroNeumannBoundary::applyLoX(Field& f)
{
  double x_lo = f.getLo(0);
  double x_lo_inner = f.getInnerLo()[0];

  double y_lo = f.getLo(1);
  double y_hi = f.getHi(1);
  for (int y=y_lo; y<=y_hi; ++y)
  {
    double val = f(x_lo_inner,y);
    for (int x=x_lo; x<x_lo_inner; ++x)
      f(x,y) = val;
  }
}

void ZeroNeumannBoundary::applyLoY(Field& f)
{
  double y_lo = f.getLo(1);
  double y_lo_inner = f.getInnerLo()[1];

  double x_lo = f.getLo(0);
  double x_hi = f.getHi(0);
  for (int x=x_lo; x<=x_hi; ++x)
  {
    double val = f(x,y_lo_inner);
    for (int y=y_lo; y<y_lo_inner; ++y)
      f(x,y) = val;
  }
}

void ZeroNeumannBoundary::applyHiX(Field& f)
{
  double x_hi = f.getHi(0);
  double x_hi_inner = f.getInnerHi()[0];

  double y_lo = f.getLo(1);
  double y_hi = f.getHi(1);
  for (int y=y_lo; y<=y_hi; ++y)
  {
    double val = f(x_hi_inner,y);
    for (int x=x_hi; x>x_hi_inner; --x)
      f(x,y) = val;
  }
}

void ZeroNeumannBoundary::applyHiY(Field& f)
{
  double y_hi = f.getHi(1);
  double y_hi_inner = f.getInnerHi()[1];

  double x_lo = f.getLo(0);
  double x_hi = f.getHi(0);
  for (int x=x_lo; x<=x_hi; ++x)
  {
    double val = f(x,y_hi_inner);
    for (int y=y_hi; y>y_hi_inner; --y)
      f(x,y) = val;
  }
}
