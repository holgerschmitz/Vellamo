/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "vellamo.hpp"
#include "diagnostic.hpp"
#include "solver.hpp"

#include <schnek/parser.hpp>
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/optional.hpp>

#include <mpi.h>

#include <fstream>
#include <string>
#include <unistd.h>


double step(double x, double x0)
{
  return (x>=x0)?1.0:0.0;
}

double stepi(double x, double x0)
{
  return (x>=x0)?0.0:1.0;
}

double box(double x, double xmin, double xmax)
{
  return ((x>=xmin)&&(x<xmax))?1.0:0.0;
}

Vellamo *Vellamo::instance;

Vellamo::Vellamo()
{
  instance = this;
}

void Vellamo::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.5);
  x_parameters = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);

  Rho_parameter = parameters.addParameter("Rho", &initRho, 0.0);

  M_parameters = parameters.addArrayParameter("M", initM, 0.0);

  E_parameter = parameters.addParameter("E", &initE, 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);
}

void Vellamo::registerData()
{
  addData("Rho", Rho);

  addData("Mx", M[0]);
  addData("My", M[1]);

  addData("E", E);
}


void Vellamo::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  updater.addIndependentArray(x_parameters);
  schnek::fill_field(*Rho, x, initRho, updater, Rho_parameter);

  schnek::fill_field(*M[0], x, initM[0], updater, M_parameters[0]);
  schnek::fill_field(*M[1], x, initM[1], updater, M_parameters[1]);

  schnek::fill_field(*E, x, initE, updater, E_parameter);
}

void Vellamo::init()
{
  globalMax = gridSize-1;

  subdivision.init(gridSize, 2);

  dx[0] = size[0] / gridSize[0];
  dx[1] = size[1] / gridSize[1];

  dt = cflFactor*std::min(dx[0],dx[1])/clight;

  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  innerRange = Range(lowIn, highIn);
  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>(0,0), size);
  schnek::Array<bool, DIMENSION> stagger;

  stagger = false;

  Rho = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  M[0] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);
  M[1] = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  E = boost::make_shared<Field>(lowIn, highIn, domainSize, stagger, 2);

  fillValues();
}

void Vellamo::execute()
{

  schnek::DiagnosticManager::instance().setTimeCounter(&timestep);

  if (numChildren()==0)
    throw schnek::VariableNotFoundException("At least one Fluid Solver needs to be specified");

  double time = 0.0;
  timestep = 0;

  while (time<=tMax)
  {
    schnek::DiagnosticManager::instance().execute();

    if (subdivision.master())
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;

    boost::optional<double> maxDt;
    BOOST_FOREACH(Solver *f, childBlocks())
    {
      maxDt = (maxDt)?std::min(maxDt.get(), f->maxDt()):f->maxDt();
    }

    dt = cflFactor*maxDt.get();

    BOOST_FOREACH(Solver *f, childBlocks())
    {
      f->rungeKuttaStep(dt);
    }

    time += dt;
    ++timestep;
  }

  schnek::DiagnosticManager::instance().execute();
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try
  {
    schnek::BlockClasses blocks;

    blocks.registerBlock("vellamo").setClass<Vellamo>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();
    blocks("CompressibleEuler").setClass<Solver>();

    blocks("vellamo").addChildren("FieldDiag")("CompressibleEuler");

    std::ifstream in("vellamo.setup");
    if (!in) throw std::string("Could not open file 'vellamo.setup'");

    schnek::Parser P("vellamo", "vellamo", blocks);
    registerCMath(P.getFunctionRegistry());

    P.getFunctionRegistry().registerFunction("step",step);
    P.getFunctionRegistry().registerFunction("stepi",stepi);
    P.getFunctionRegistry().registerFunction("box",box);

    schnek::pBlock application = P.parse(in);

    Vellamo &mpulse = dynamic_cast<Vellamo&>(*application);
    mpulse.initAll();

    if (mpulse.getSubdivision().master())
    {
      std::ofstream referencesText("information.tex");
      std::ofstream referencesBib("references.bib");

      schnek::LiteratureManager::instance().writeInformation(referencesText,"references.bib");
      schnek::LiteratureManager::instance().writeBibTex(referencesBib);
      referencesText.close();
      referencesBib.close();
    }

    mpulse.execute();
  }
  catch (schnek::ParserError &e)
  {
    std::cerr << "Parse error in " << e.getFilename() << " at line "
        << e.getLine() << ": " << e.message << std::endl;
    return -1;
  }
  catch (schnek::VariableNotInitialisedException &e)
  {
    std::cerr << "Variable was not initialised: " << e.getVarName() << std::endl;
    return -1;
  }
  catch (schnek::EvaluationException &e)
  {
    std::cerr << "Error in evaluation: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (schnek::VariableNotFoundException &e)
  {
    std::cerr << "Error: " << e.getMessage() << std::endl;
    return -1;
  }
  catch (SchnekException &e)
  {
    std::cerr << "An error occured" << std::endl;
    return -1;
  }
  catch (std::string &err)
  {
    std::cerr << "FATAL ERROR: " << err << std::endl;
    return -1;
  }

  MPI_Finalize();

  return 0;
}

// end of main
