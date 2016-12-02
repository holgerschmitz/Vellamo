/*
 * solver.cpp
 *
 *  Created on: 7 Oct 2016
 *      Author: Holger Schmitz
 */

#include "euler_solver.hpp"

#include "boundary.hpp"

#include <schnek/tools/literature.hpp>

#include <boost/make_shared.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>

static const double adiabaticGamma = 1.4; //5.0/3.0;

void AdiabaticSolver::init()
{
  retrieveData("Rho", Rho);

  retrieveData("Mx", Mx);
  retrieveData("My", My);

  retrieveData("E", E);

  schnek::LiteratureArticle Kurganov2000("Kurganov2000", "A. Kurganov and S. Noelle and G. Petrova",
      "Semidiscrete central-upwind schemes for hyperbolic conservation laws and Hamilton--Jacobi equations",
      "SIAM J. Sci. Comput.", "2001", "23", "707");

  schnek::LiteratureManager::instance().addReference(
      "Semidiscrete central-upwind scheme for hyperbolic conservation laws", Kurganov2000);
}

void AdiabaticSolver::postInit()
{
  Rho_s = boost::make_shared<Field>(*Rho);
  Mx_s =  boost::make_shared<Field>(*Mx);
  My_s =  boost::make_shared<Field>(*My);
  E_s  =  boost::make_shared<Field>(*E);

  dx = Vellamo::getDx();


}

inline void AdiabaticSolver::checkFluid(const FluidValues &u)
{
  return;

  for (int i=0; i<4; ++i)
    if (!(u[i]==u[i]) || !(0.0*u[i]==0.0*u[i]))
    {
      std::cerr << "NaN or Infinity\n";
    }
}

inline double AdiabaticSolver::van_leer(double u, double up, double um)
{
  double du = (up-u)*(u-um);

  return (du>0.0)?du/(up-um):0.0;
}

inline double AdiabaticSolver::speed_cf(double rho, double p)
{
 return (p>0.0)?(0.5*sqrt(4.0*adiabaticGamma*p/rho)):0.0;
}


inline double AdiabaticSolver::eqn_state_ideal_gas(FluidValues& u)
{
  double internal_energy = std::max(0.0, u[C_E] - 0.5*(u[C_MX]*u[C_MX] + u[C_MY]*u[C_MY])/u[C_RHO]);

  return (adiabaticGamma-1.0)*internal_energy;
}

void AdiabaticSolver::minmax_local_speed(int d, FluidValues uW, FluidValues uE, double pW, double pE, double &ap, double &am)
{
  double vW, vE;
  double cfW, cfE;

  vW = uW[C_MX+d] / uW[C_RHO];
  vE = uE[C_MX+d] / uE[C_RHO];

  cfW = speed_cf(uW[C_RHO], pW);
  cfE = speed_cf(uE[C_RHO], pE);

  ap = std::max( (vW+cfW), std::max( (vE+cfE), 0.0 ));
  am = std::min( (vW-cfW), std::min( (vE-cfE), 0.0 ));
}

void AdiabaticSolver::reconstruct_x(int i, int j, int dir, FluidValues& u)
{
  u[C_RHO] = (*Rho)(i,j) + dir*van_leer((*Rho)(i,j), (*Rho)(i+1,j), (*Rho)(i-1,j));
  u[C_MX]  = (*Mx)(i,j)  + dir*van_leer((*Mx)(i,j),  (*Mx)(i+1,j),  (*Mx)(i-1,j));
  u[C_MY]  = (*My)(i,j)  + dir*van_leer((*My)(i,j),  (*My)(i+1,j),  (*My)(i-1,j));
  u[C_E]   = (*E)(i,j)   + dir*van_leer((*E)(i,j),   (*E)(i+1,j),   (*E)(i-1,j));
}

void AdiabaticSolver::reconstruct_y(int i, int j, int dir, FluidValues& u)
{
  u[C_RHO] = (*Rho)(i,j) + dir*van_leer((*Rho)(i,j), (*Rho)(i,j+1), (*Rho)(i,j-1));
  u[C_MX]  = (*Mx)(i,j)  + dir*van_leer((*Mx)(i,j),  (*Mx)(i,j+1),  (*Mx)(i,j-1));
  u[C_MY]  = (*My)(i,j)  + dir*van_leer((*My)(i,j),  (*My)(i,j+1),  (*My)(i,j-1));
  u[C_E]   = (*E)(i,j)   + dir*van_leer((*E)(i,j),   (*E)(i,j+1),   (*E)(i,j-1));
}

void AdiabaticSolver::flux_function_x(FluidValues u, double p, FluidValues &f)
{
  double rho = u[C_RHO];
  double mx = u[C_MX];
  double my = u[C_MY];
  double engy = u[C_E];

  f[C_RHO]   = mx;
  f[C_E]     = (engy + p)*mx/rho;
  f[C_MX]    = mx*mx/rho + p;
  f[C_MY]    = mx*my/rho;
}

void AdiabaticSolver::flux_function_y(FluidValues u, double p, FluidValues &f)
{
  double rho = u[C_RHO];
  double mx = u[C_MX];
  double my = u[C_MY];
  double engy = u[C_E];

  f[C_RHO]   = my;
  f[C_E]     = (engy + p)*my/rho;
  f[C_MX]    = mx*my/rho;
  f[C_MY]    = my*my/rho + p;
}

inline void AdiabaticSolver::flux_x(int i, int j, FluidValues& flux)
{
  FluidValues uW, uE;
  double ap, am;
  FluidValues fE, fW;
  double pE, pW;

  // reconstruct the MHD variables on the cell boundary
  reconstruct_x(i,   j, +1, uE);
  reconstruct_x(i+1, j, -1, uW);

  // calculate the thermodynamical variables, pressure and temperature
  pE = eqn_state_ideal_gas(uE);
  pW = eqn_state_ideal_gas(uW);

  // determine the minimum and maximum local speeds
  minmax_local_speed(0, uW, uE, pW, pE, ap, am);

//  if (ap==am) {
//    flux = 0.0;
//    return;
//  }

  // evaluate the flux function
  flux_function_x(uW, pW, fW);
  flux_function_x(uE, pE, fE);

  // assemble everything to calculate the flux
  flux = (ap*fE - am*fW + ap*am*(uW-uE))/(ap-am);
  checkFluid(flux);
}

inline void AdiabaticSolver::flux_y(int i, int j, FluidValues& flux)
{
  FluidValues uW, uE;
  double ap, am;
  FluidValues fE, fW;
  double pE, pW;

  // reconstruct the MHD variables on the cell boundary
  reconstruct_y(i, j  , +1, uE);
  reconstruct_y(i, j+1, -1, uW);

  // calculate the thermodynamical variables, pressure and temperature
  pE = eqn_state_ideal_gas(uE);
  pW = eqn_state_ideal_gas(uW);

  // determine the minimum and maximum local speeds
  minmax_local_speed(1, uW, uE, pW, pE, ap, am);

//  if (ap==am) {
//    flux = 0.0;
//    return;
//  }

  // evaluate the flux function
  flux_function_y(uW, pW, fW);
  flux_function_y(uE, pE, fE);

  // assemble everything to calculate the flux
  flux = (ap*fE-am*fW + ap*am*(uW-uE))/(ap-am);
  checkFluid(flux);
}

inline void AdiabaticSolver::hydroRhs(Index p, FluidValues& dudt)
{
  FluidValues fmx, fpx;
  FluidValues fmy, fpy;

  flux_x(p[0]-1, p[1], fmx);
  flux_x(  p[0], p[1], fpx);

  flux_y(p[0], p[1]-1, fmy);
  flux_y(p[0],   p[1], fpy);

  dudt = (fmx-fpx)/dx[0] + (fmy-fpy)/dx[1];
  checkFluid(dudt);
}


void AdiabaticSolver::timeStep(double dt)
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  Field &Rho = *(this->Rho);
  Field &Mx = *(this->Mx);
  Field &My = *(this->My);
  Field &E = *(this->E);

  Field &Rho_s = *(this->Rho_s);
  Field &Mx_s = *(this->Mx_s);
  Field &My_s = *(this->My_s);
  Field &E_s = *(this->E_s);

  Index lo = Rho.getInnerLo();
  Index hi = Rho.getInnerHi();
  Index p;
  FluidValues dudt;

  // First step
  for (p[0]=lo[0]; p[0]<=hi[0]; ++p[0])
    for (p[1]=lo[1]; p[1]<=hi[1]; ++p[1])
    {
      hydroRhs(p, dudt);

      Rho_s(p[0],p[1]) = Rho(p[0],p[1]) + dt*dudt[C_RHO];
      Mx_s(p[0],p[1])  = Mx(p[0],p[1])  + dt*dudt[C_MX];
      My_s(p[0],p[1])  = My(p[0],p[1])  + dt*dudt[C_MY];
      E_s(p[0],p[1])   = E(p[0],p[1])   + dt*dudt[C_E];
    }

  // Swap starred fields and the unstarred arrays
  for (p[0]=lo[0]; p[0]<=hi[0]; ++p[0])
      for (p[1]=lo[1]; p[1]<=hi[1]; ++p[1])
      {
        std::swap(Rho(p[0],p[1]), Rho_s(p[0],p[1]));
        std::swap(Mx(p[0],p[1]) , Mx_s(p[0],p[1]));
        std::swap(My(p[0],p[1]) , My_s(p[0],p[1]));
        std::swap(E(p[0],p[1])  , E_s(p[0],p[1]));
      }

  subdivision.exchange(Rho);
  subdivision.exchange(Mx);
  subdivision.exchange(My);
  subdivision.exchange(E);

  BOOST_FOREACH(pBoundaryCondition bc, schnek::BlockContainer<BoundaryCondition>::childBlocks())
  {
    bc->apply(Rho, Mx, My, E);
  }

  // second step
  for (p[0]=lo[0]; p[0]<=hi[0]; ++p[0])
    for (p[1]=lo[1]; p[1]<=hi[1]; ++p[1])
    {
      hydroRhs(p, dudt);

      Rho_s(p[0],p[1]) = 0.5*(Rho(p[0],p[1]) + Rho_s(p[0],p[1]) + dt*dudt[C_RHO]);
      Mx_s(p[0],p[1])  = 0.5*(Mx(p[0],p[1])  + Mx_s(p[0],p[1])  + dt*dudt[C_MX]);
      My_s(p[0],p[1])  = 0.5*(My(p[0],p[1])  + My_s(p[0],p[1])  + dt*dudt[C_MY]);
      E_s(p[0],p[1])   = 0.5*(E(p[0],p[1])   + E_s(p[0],p[1])   + dt*dudt[C_E]);
    }


  // Copy starred back into unstarred fields

  for (p[0]=lo[0]; p[0]<=hi[0]; ++p[0])
    for (p[1]=lo[1]; p[1]<=hi[1]; ++p[1])
    {
      Rho(p[0],p[1]) = Rho_s(p[0],p[1]);
      Mx(p[0],p[1]) = Mx_s(p[0],p[1]);
      My(p[0],p[1]) = My_s(p[0],p[1]);
      E(p[0],p[1]) = E_s(p[0],p[1]);
    }

  subdivision.exchange(Rho);
  subdivision.exchange(Mx);
  subdivision.exchange(My);
  subdivision.exchange(E);

  BOOST_FOREACH(pBoundaryCondition bc, schnek::BlockContainer<BoundaryCondition>::childBlocks())
  {
    bc->apply(Rho, Mx, My, E);
  }
}

double AdiabaticSolver::maxDt()
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  Field &Rho = *(this->Rho);
  Field &Mx = *(this->Mx);
  Field &My = *(this->My);
  Field &E = *(this->E);

  Index lo = Rho.getInnerLo();
  Index hi = Rho.getInnerHi();
  FluidValues u;
  double p;
  double max_speed = 0.0;

  double min_dx = std::min(dx[0], dx[1]);

  for (int i=lo[0]; i<=hi[0]; ++i)
    for (int j=lo[1]; j<=hi[1]; ++j)
    {
      u[C_RHO]    = Rho(i,j);
      u[C_MX]    = Mx(i,j);
      u[C_MY]    = My(i,j);
      u[C_E] = E(i,j);

      p = eqn_state_ideal_gas(u);

      double v_max = std::max(fabs(u[C_MX]),fabs(u[C_MY]))/u[C_RHO];

      max_speed = std::max(max_speed,speed_cf(u[C_RHO], p)+v_max);
    }

  max_speed = subdivision.maxReduce(max_speed);
  return min_dx/max_speed;
}
