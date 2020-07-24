/*
 * fieldsolver.hpp
 *
 *  Created on: 5 Feb 2008
 *      Author: Holger Schmitz
 */

#ifndef EULER_SOLVER_HPP_
#define EULER_SOLVER_HPP_

#include "vellamo.hpp"
//#include "solver.hpp"
//
//#include <schnek/variables/blockcontainer.hpp>
//
//class Storage;
//class Current;
//
//const int C_RHO = 0;
//const int C_E   = 1;
//const int C_M[]  = {2, 3, 4};
//
//class AdiabaticSolver: public Solver
//{
//  public:
//    typedef schnek::Array<double, DIMENSION+2> FluidValues;
//  private:
//    pField Rho;
//    schnek::Array<pField, DIMENSION> M;
//    pField E;
//
//    pField Rho_s;
//    schnek::Array<pField, DIMENSION> M_s;
//    pField E_s;
//
//    Vector dx;
//
//    double van_leer(double u, double up, double um);
//    double speed_cf(double rho, double p);
//    double eqn_state_ideal_gas(FluidValues &u);
//    void minmax_local_speed(size_t dim, FluidValues uW, FluidValues uE, double pW, double pE, double &ap, double &am);
//    void reconstruct(size_t dim, const Index &pos, int dir, FluidValues &u);
//    void reconstruct_y(int i, int j, int dir, FluidValues &u);
//    void flux_function(size_t dim, FluidValues u, double p, FluidValues &f);
//    void flux_function_y(FluidValues u, double p, FluidValues &f);
//    void flux(size_t dim, const Index &pos, FluidValues& flux);
//    void flux_y(int i, int j, FluidValues &flux);
//    void hydroRhs(Index p, FluidValues &dudt);
//
//    void checkFluid(const FluidValues &dudt);
//  public:
//    void init();
//    void postInit();
//    void timeStep(double dt);
//    double maxDt();
//};

#endif // EULER_SOLVER_HPP_
