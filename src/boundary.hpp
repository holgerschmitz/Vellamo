/*
 * boundary.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "vellamo.hpp"


//class ZeroNeumannBoundary
//{
//  public:
//    void applyLo(size_t dim, Field &f);
//    void applyHi(size_t dim, Field &f);
//};
//
//class ZeroDirichletBoundary
//{
//  public:
//    void applyLo(size_t dim, Field &f);
//    void applyHi(size_t dim, Field &f);
//};
//
//class BoundaryCondition : public schnek::ChildBlock<BoundaryCondition>
//{
//  private:
//    Index applyLo;
//    Index applyHi;
//  public:
//    void initParameters(schnek::BlockParameters &blockPars);
//    virtual ~BoundaryCondition() {}
//    void apply(Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E);
//
//    virtual void applyLoDim(size_t size_t, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E) = 0;
//    virtual void applyHiDim(size_t size_t, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E) = 0;
//};
//
//typedef boost::shared_ptr<BoundaryCondition> pBoundaryCondition;
//
//class ZeroNeumannBoundaryBlock : public BoundaryCondition
//{
//  private:
//    ZeroNeumannBoundary boundary;
//  public:
//    void applyLoDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E);
//    void applyHiDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E);
//};
//
//class WallBoundaryBlock : public BoundaryCondition
//{
//  private:
//    ZeroNeumannBoundary neumannBoundary;
//    ZeroDirichletBoundary dirichletBoundary;
//  public:
//    void applyLoDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E);
//    void applyHiDim(size_t dim, Field &Rho, schnek::Array<pField, DIMENSION> M, Field &E);
//};


#endif /* BOUNDARY_HPP_ */
