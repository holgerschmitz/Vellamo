/*
 * vellamo.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef VELLAMO_HPP_
#define VELLAMO_HPP_

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "types.hpp"

class Solver;
class HydroFields;

class Vellamo : public schnek::Block,
                public schnek::BlockContainer<Solver>,
                public schnek::BlockContainer<HydroFields>,
                public boost::enable_shared_from_this<Vellamo>
{
  private:
    static Vellamo *instance;
    Index globalMax;
    Index gridSize;
    Vector size;
    Vector dx;

    double cflFactor;
    double dt;
    schnek::MPICartSubdivision<Field> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;
    schnek::pParametersGroup spaceVars;

    double tMax;
    Range innerRange;

    int timestep;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void initFields();
  public:
    Vellamo();
    void init();
    void execute();

    static Index getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static Vector getSize() { return instance->size; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
    static Vector &getX() { return instance->x; }
    static schnek::Array<schnek::pParameter, DIMENSION> &getXParameter() { return instance->x_parameters; }
};

#endif // VELLAMO_HPP_
