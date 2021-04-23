/*
 * vellamo.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef VELLAMO_HPP_
#define VELLAMO_HPP_

#include "../huerto/types.hpp"
#include "../huerto/constants.hpp"
#include "../huerto/hydrodynamics/hydro_fields.hpp"
#include "../huerto/hydrodynamics/hydro_solver.hpp"
#include "../huerto/hydrodynamics/euler/euler_knp.hpp"
#include "../huerto/hydrodynamics/euler/adiabatic_knp.hpp"

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

class Vellamo : public schnek::Block,
                public schnek::BlockContainer<HydroSolver>,
                public schnek::BlockContainer<HydroFields>,
                public boost::enable_shared_from_this<Vellamo>,
                public SimulationContext
{
  private:
    double cflFactor;

    schnek::pParametersGroup spaceVars;

    Range innerRange;

    int timestep;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void initFields();
  public:
    Vellamo();
    void init();
    void execute();

    schnek::Array<schnek::pParameter, DIMENSION> &getXParameter() { return x_parameters; }
};

#endif // VELLAMO_HPP_
