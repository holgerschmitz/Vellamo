/*
 * fieldsolver.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_EULER_SOLVER_H
#define MPULSE_EULER_SOLVER_H

#include "vellamo.hpp"
#include "solver.hpp"

#include <schnek/variables/blockcontainer.hpp>

class Storage;
class Current;

const int C_RHO = 0;
const int C_MX  = 1;
const int C_MY  = 2;
const int C_E   = 3;

class EulerSolver: public Solver
{
  public:
    typedef schnek::Array<double, 4> FluidValues;
  private:
    pField Rho;
    pField Mx, My;
    pField E;

    pField Rho_s;
    pField Mx_s, My_s;
    pField E_s;

    Vector dx;

    double van_leer(double u, double up, double um);
    double speed_cf(double rho, double p);
    double eqn_state_ideal_gas(FluidValues &u);
    void minmax_local_speed(int d, FluidValues uW, FluidValues uE, double pW, double pE, double &ap, double &am);
    void reconstruct_x(int i, int j, int dir, FluidValues &u);
    void reconstruct_y(int i, int j, int dir, FluidValues &u);
    void flux_function_x(FluidValues u, double p, FluidValues &f);
    void flux_function_y(FluidValues u, double p, FluidValues &f);
    void flux_x(int i, int j, FluidValues &flux);
    void flux_y(int i, int j, FluidValues &flux);
    void hydroRhs(Index p, FluidValues &dudt);

    void checkFluid(const FluidValues &dudt);
  public:
    void init();
    void postInit();
    void timeStep(double dt);
    double maxDt();
};

#endif // MPULSE_EULER_SOLVER_H
