/*
 * boundary.cpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#include "boundary.hpp"


void ZeroNeumannBoundary::applyLo(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerLo = f.getInnerLo()[dim];
  range.getHi()[dim] =  innerLo - 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    Index src = *it;
    src[dim] = innerLo;
    f[*it] = f[src];
  }
}


void ZeroNeumannBoundary::applyHi(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerHi = f.getInnerHi()[dim];
  range.getLo()[dim] =  innerHi + 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    Index src = *it;
    src[dim] = innerHi;
    f[*it] = f[src];
  }
}


void ZeroDirichletBoundary::applyLo(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerLo = f.getInnerLo()[dim];
  range.getHi()[dim] =  innerLo - 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    f[*it] = 0.0;
  }
}

void ZeroDirichletBoundary::applyHi(size_t dim, Field& f)
{
  Index lo_inner = f.getInnerLo();

  Range range(f.getLo(), f.getHi());

  int innerHi = f.getInnerHi()[dim];
  range.getLo()[dim] =  innerHi + 1;

  const Range::iterator itEnd = range.end();

  for (Range::iterator it = range.begin();
       it != itEnd;
       ++it)
  {
    f[*it] = 0.0;
  }
}

void BoundaryCondition::initParameters(schnek::BlockParameters &blockPars)
{
  blockPars.addArrayParameter("low_",applyLo, Index(0));
  blockPars.addArrayParameter("high_",applyHi, Index(0));
}

void BoundaryCondition::apply(Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E)
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  for (size_t i=0; i<DIMENSION; i++)
  {
    if (bool(applyLo[i]) && subdivision.isBoundLo(i)) applyLoDim(i, Rho, M, E);
    if (bool(applyHi[i]) && subdivision.isBoundHi(i)) applyHiDim(i, Rho, M, E);
  }
}

void ZeroNeumannBoundaryBlock::applyLoDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E)
{
  boundary.applyLo(dim, Rho);
  for (size_t i=0; i<DIMENSION; i++)
  {
    boundary.applyLo(dim, *M[i]);
  }
  boundary.applyLo(dim, E);
}


void ZeroNeumannBoundaryBlock::applyHiDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E)
{
  boundary.applyHi(dim, Rho);
  for (size_t i=0; i<DIMENSION; i++)
  {
    boundary.applyHi(dim, *M[i]);
  }
  boundary.applyHi(dim, E);
}

void WallBoundaryBlock::applyLoDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E)
{
  neumannBoundary.applyLo(dim, Rho);
  for (size_t i=0; i<DIMENSION; i++)
  {
    if (i == dim) {
      dirichletBoundary.applyLo(dim, *M[i]);
    } else {
      neumannBoundary.applyLo(dim, *M[i]);
    }
  }
  neumannBoundary.applyLo(dim, E);
}


void WallBoundaryBlock::applyHiDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E)
{
  neumannBoundary.applyHi(dim, Rho);
  for (size_t i=0; i<DIMENSION; i++)
  {
    if (i == dim)
    {
      dirichletBoundary.applyHi(dim, *M[i]);
    } else {
      neumannBoundary.applyHi(dim, *M[i]);
    }
  }
  neumannBoundary.applyHi(dim, E);
}
