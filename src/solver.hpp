/*
 * solver.hpp
 *
 *  Created on: 14 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

class Solver : public schnek::ChildBlock<Solver>
{
  public:
    virtual ~Solver() {}
    virtual double maxDt() = 0;
    virtual void timeStep(double dt) = 0;
};



#endif /* SOLVER_HPP_ */
