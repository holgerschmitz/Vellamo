/*
 * em_fields.cpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */

#include "hydro_fields.hpp"
#include "vellamo.hpp"

#include <schnek/grid/domainsubdivision.hpp>
#include <schnek/tools/fieldtools.hpp>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include <string>
#include <iostream>

void HydroFields::registerData()
{
  Rho = boost::make_shared<Field>();
  addData("Rho", Rho);

  const std::string coords[] = {"x", "y", "z"};
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i] = boost::make_shared<Field>();
    addData(std::string("M")+coords[i], M[i]);
  }

  E = boost::make_shared<Field>();
  addData("E", E);
}

void HydroFields::initParameters(schnek::BlockParameters &parameters)
{
  Rho_parameter = parameters.addParameter("Rho", &initRho, 0.0);
  M_parameters = parameters.addArrayParameter("M", initM, 0.0);
  E_parameter = parameters.addParameter("E", &initE, 0.0);
}

void HydroFields::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  Vector &x = Vellamo::getX();
  schnek::Array<schnek::pParameter, DIMENSION> x_parameters = Vellamo::getXParameter();

  updater.addIndependentArray(x_parameters);
  schnek::fill_field(*Rho, x, initRho, updater, Rho_parameter);

  schnek::fill_field(*M[0], x, initM[0], updater, M_parameters[0]);
  schnek::fill_field(*M[1], x, initM[1], updater, M_parameters[1]);

  schnek::fill_field(*E, x, initE, updater, E_parameter);
}

void HydroFields::init()
{
  const schnek::DomainSubdivision<Field> &subdivision = Vellamo::getSubdivision();
  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>::Zero(), Vellamo::getSize());

  schnek::Array<bool, DIMENSION> stagger;
  stagger = false;

  Rho->resize(lowIn, highIn, domainSize, stagger, 2);
  for (size_t i=0; i<DIMENSION; ++i)
  {
    M[i]->resize(lowIn, highIn, domainSize, stagger, 2);
  }
  E->resize(lowIn, highIn, domainSize, stagger, 2);
}
