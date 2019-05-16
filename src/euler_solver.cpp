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

  const std::string coords[] = {"x", "y", "z"};
  for (size_t i=0; i<DIMENSION; ++i)
  {
    retrieveData(std::string("M")+coords[i], M[i]);
  }

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
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M_s[i] =  boost::make_shared<Field>(*M[i]);
  }

  dx = Vellamo::getDx();


}

inline void AdiabaticSolver::checkFluid(const FluidValues &u)
{
  return;

  for (size_t i=0; i<4; ++i)
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
  double sqrU = 0.0;
  for (size_t i=0; i<DIMENSION; ++i)
  {
    sqrU += u[C_M[i]]*u[C_M[i]];
  }

  double internal_energy = std::max(0.0, u[C_E] - 0.5*(sqrU)/u[C_RHO]);

  return (adiabaticGamma-1.0)*internal_energy;
}

void AdiabaticSolver::minmax_local_speed(size_t dim, FluidValues uW, FluidValues uE, double pW, double pE, double &ap, double &am)
{
  double vW, vE;
  double cfW, cfE;

  vW = uW[C_M[dim]] / uW[C_RHO];
  vE = uE[C_M[dim]] / uE[C_RHO];

  cfW = speed_cf(uW[C_RHO], pW);
  cfE = speed_cf(uE[C_RHO], pE);

  ap = std::max( (vW+cfW), std::max( (vE+cfE), 0.0 ));
  am = std::min( (vW-cfW), std::min( (vE-cfE), 0.0 ));
}

void AdiabaticSolver::reconstruct(size_t dim, const Index &pos, int dir, FluidValues& u)
{
  Index posp = pos;
  ++posp[dim];
  Index posm = pos;
  --posm[dim];

  u[C_RHO] = (*Rho)[pos] + dir*van_leer((*Rho)[pos], (*Rho)[posp], (*Rho)[posm]);

  for (size_t i=0; i<DIMENSION; ++i)
  {
    u[C_M[i]] = (*M[i])[pos] + dir*van_leer((*M[i])[pos],  (*M[i])[posp],  (*M[i])[posm]);
  }
  u[C_E] = (*E)[pos] + dir*van_leer((*E)[pos],   (*E)[posp],   (*E)[posm]);
}

void AdiabaticSolver::flux_function(size_t dim, FluidValues u, double p, FluidValues &f)
{
  double rho = u[C_RHO];
  double mdim = u[C_M[dim]];
  double engy = u[C_E];

  f[C_RHO]   = mdim;
  f[C_E]     = (engy + p)*mdim/rho;
  for (size_t i=0; i<DIMENSION; ++i)
  {
    f[C_M[i]] = mdim*u[C_M[i]]/rho;
  }
  f[C_M[dim]] += p;
}


inline void AdiabaticSolver::flux(size_t dim, const Index &pos, FluidValues& flux)
{
  FluidValues uW, uE;
  double ap, am;
  FluidValues fE, fW;
  double pE, pW;

  Index posp = pos;
  ++posp[dim];

  // reconstruct the MHD variables on the cell boundary
  reconstruct(dim, pos,  +1, uE);
  reconstruct(dim, posp, -1, uW);

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
  flux_function(dim, uW, pW, fW);
  flux_function(dim, uE, pE, fE);

  // assemble everything to calculate the flux
  flux = (ap*fE - am*fW + ap*am*(uW-uE))/(ap-am);
  checkFluid(flux);
}

inline void AdiabaticSolver::hydroRhs(Index pos, FluidValues& dudt)
{

  FluidValues sum = 0;
  for (size_t i=0; i<DIMENSION; ++i)
  {
    Index posm = pos;
    FluidValues fm, fp;
    --posm[i];
    flux(i, posm, fm);
    flux(i, pos,  fp);

    sum += (fm - fp) / dx[i];
  }

  dudt = sum;
  checkFluid(dudt);
}


void AdiabaticSolver::timeStep(double dt)
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  Field &Rho = *(this->Rho);
  schnek::Array<Field*, DIMENSION> M;

  Field &E = *(this->E);

  Field &Rho_s = *(this->Rho_s);
  schnek::Array<Field*, DIMENSION> M_s;
  Field &E_s = *(this->E_s);

  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i] = &(*this->M[i]);
    M_s[i] = &(*this->M_s[i]);
  }

  Index lo = Rho.getInnerLo();
  Index hi = Rho.getInnerHi();
  Range range(lo, hi);
  Range::iterator range_end = range.end();
  FluidValues dudt;

  // First step
  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    hydroRhs(p, dudt);

    Rho_s[p] = Rho[p] + dt*dudt[C_RHO];
    for (size_t i=0; i<DIMENSION; ++i)
    {
     (*M_s[i])[p] = (*M[i])[p]  + dt*dudt[C_M[i]];
    }
    E_s[p]   = E[p]   + dt*dudt[C_E];
  }

  // Swap starred fields and the unstarred arrays
  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    std::swap(Rho[p], Rho_s[p]);
    for (size_t i=0; i<DIMENSION; ++i)
    {
      std::swap((*M[i])[p], (*M_s[i])[p]);
    }
    std::swap(E[p], E_s[p]);
  }

  subdivision.exchange(Rho);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    subdivision.exchange(*M[i]);
  }
  subdivision.exchange(E);

  BOOST_FOREACH(pBoundaryCondition bc, schnek::BlockContainer<BoundaryCondition>::childBlocks())
  {
    bc->apply(Rho, this->M, E);
  }

  // second step
  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    hydroRhs(p, dudt);

    Rho_s[p] = 0.5*(Rho[p] + Rho_s[p] + dt*dudt[C_RHO]);
    for (size_t i=0; i<DIMENSION; ++i)
    {
     (*M_s[i])[p] = 0.5*((*M[i])[p]  + (*M_s[i])[p]  + dt*dudt[C_M[i]]);
    }
    E_s[p]   = 0.5*(E[p]   + E_s[p]   + dt*dudt[C_E]);
  }


  // Copy starred back into unstarred fields

  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    Rho[p] = Rho_s[p];
    for (size_t i=0; i<DIMENSION; ++i)
    {
     (*M[i])[p] = (*M_s[i])[p];
    }
    E[p] = E_s[p];
  }

  subdivision.exchange(Rho);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    subdivision.exchange(*M[i]);
  }
  subdivision.exchange(E);

  BOOST_FOREACH(pBoundaryCondition bc, schnek::BlockContainer<BoundaryCondition>::childBlocks())
  {
    bc->apply(Rho, this->M, E);
  }
}

double AdiabaticSolver::maxDt()
{
  schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();

  Field &Rho = *(this->Rho);

  schnek::Array<Field*, DIMENSION> M;
  Field &E = *(this->E);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i] = &(*this->M[i]);
  }

  Index lo = Rho.getInnerLo();
  Index hi = Rho.getInnerHi();
  Range range(lo, hi);
  Range::iterator range_end = range.end();

  FluidValues u;
  double max_speed = 0.0;

  double min_dx = std::min(dx[0], dx[1]);

  for (Range::iterator it = range.begin();
       it != range_end;
       ++it)
  {
    const Index &p = *it;
    u[C_RHO]    = Rho[p];

    double maxU = 0.0;
    for (size_t i=0; i<DIMENSION; ++i)
    {
      u[C_M[i]]    = (*M[i])[p];
      maxU = std::max(maxU, fabs(u[C_M[i]]));
    }
    u[C_E] = E[p];

    double pressure = eqn_state_ideal_gas(u);

    double v_max = maxU/u[C_RHO];

    max_speed = std::max(max_speed,speed_cf(u[C_RHO], pressure)+v_max);
  }

  max_speed = subdivision.maxReduce(max_speed);
  return min_dx/max_speed;
}
