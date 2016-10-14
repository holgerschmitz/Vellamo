/*
 * boundary.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include "vellamo.hpp"

class BoundaryCondition : public schnek::ChildBlock<BoundaryCondition>
{
  private:
    Index applyLo;
    Index applyHi;
  public:
    void initParameters(schnek::BlockParameters &blockPars);
    virtual ~BoundaryCondition() {}
    void apply(Field &f);

    virtual void applyLoX(Field &f) = 0;
    virtual void applyLoY(Field &f) = 0;
    virtual void applyHiX(Field &f) = 0;
    virtual void applyHiY(Field &f) = 0;
};

class ZeroNeumannBoundary : public BoundaryCondition
{
  public:
    void applyLoX(Field &f);
    void applyLoY(Field &f);
    void applyHiX(Field &f);
    void applyHiY(Field &f);
};


#endif /* BOUNDARY_HPP_ */
