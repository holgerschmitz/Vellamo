/*
 * boundary.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "vellamo.hpp"


class ZeroNeumannBoundary
{
  public:
    void applyLoX(Field &f);
    void applyLoY(Field &f);
    void applyHiX(Field &f);
    void applyHiY(Field &f);
};

class ZeroDirichletBoundary
{
  public:
    void applyLoX(Field &f);
    void applyLoY(Field &f);
    void applyHiX(Field &f);
    void applyHiY(Field &f);
};

class BoundaryCondition : public schnek::ChildBlock<BoundaryCondition>
{
  private:
    Index applyLo;
    Index applyHi;
  public:
    void initParameters(schnek::BlockParameters &blockPars);
    virtual ~BoundaryCondition() {}
    void apply(Field &Rho, Field &Mx, Field &My, Field &E);

    virtual void applyLoX(Field &Rho, Field &Mx, Field &My, Field &E) = 0;
    virtual void applyLoY(Field &Rho, Field &Mx, Field &My, Field &E) = 0;
    virtual void applyHiX(Field &Rho, Field &Mx, Field &My, Field &E) = 0;
    virtual void applyHiY(Field &Rho, Field &Mx, Field &My, Field &E) = 0;
};

typedef boost::shared_ptr<BoundaryCondition> pBoundaryCondition;

class ZeroNeumannBoundaryBlock : public BoundaryCondition
{
  private:
    ZeroNeumannBoundary boundary;
  public:
    void applyLoX(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyLoY(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyHiX(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyHiY(Field &Rho, Field &Mx, Field &My, Field &E);
};

class WallBoundaryBlock : public BoundaryCondition
{
  private:
    ZeroNeumannBoundary neumannBoundary;
    ZeroDirichletBoundary dirichletBoundary;
  public:
    void applyLoX(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyLoY(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyHiX(Field &Rho, Field &Mx, Field &My, Field &E);
    void applyHiY(Field &Rho, Field &Mx, Field &My, Field &E);
};


#endif /* BOUNDARY_HPP_ */
