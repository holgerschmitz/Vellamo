/*
 * vellamo.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "vellamo.hpp"
#include "diagnostic.hpp"

#include "../huerto/hydrodynamics/euler/euler_knp.hpp"
#include "../huerto/hydrodynamics/euler/adiabatic_knp.hpp"
#include "../huerto/hydrodynamics/hydro_fields.hpp"
#include "../huerto/boundary/boundary.hpp"
#include "../huerto/constants.hpp"

#include <schnek/parser.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>
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

Vellamo::Vellamo()
{
}

void Vellamo::initParameters(schnek::BlockParameters &parameters)
{
  SimulationContext::initParameters(parameters);
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.5);
  initConstantParameters(parameters);
}


void Vellamo::initFields()
{
  if (schnek::BlockContainer<HydroFields>::childBlocks().empty())
  {
    boost::shared_ptr<HydroFields> fields(new HydroFields(shared_from_this()));
    Block::addChild(fields);
    fields->registerData();
    fields->preInit();
  }
}

void Vellamo::init()
{
  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);

  subdivision = std::make_shared<schnek::MPICartSubdivision<Field> >();
  subdivision->init(gridSize, 2);

  for (size_t i=0; i<DIMENSION; ++i)
  {
    dx[i] = size[i] / gridSize[i];
  }

  Index low  = subdivision->getLo();
  Index high = subdivision->getHi();

  Index lowIn  = subdivision->getInnerLo();
  Index highIn = subdivision->getInnerHi();

  innerRange = Range(lowIn, highIn);
  initFields();
}

void Vellamo::execute()
{
  if (schnek::BlockContainer<HydroSolver>::numChildren()==0)
    throw schnek::VariableNotFoundException("At least one Fluid Solver needs to be specified");

  time = 0.0;
  timestep = 0;
  schnek::DiagnosticManager::instance().setTimeCounter(&timestep);
  schnek::DiagnosticManager::instance().setPhysicalTime(&time);

  while (time<=tMax)
  {
    schnek::DiagnosticManager::instance().execute();

    double maxDt = std::numeric_limits<double>::max();
    for(pHydroSolver f: schnek::BlockContainer<HydroSolver>::childBlocks())
    {
      maxDt = std::min(maxDt, f->maxDt());
    }

    dt = cflFactor*maxDt;
    dt = schnek::DiagnosticManager::instance().adjustDeltaT(dt);

    if (subdivision->master())
      schnek::Logger::instance().out() <<"Time "<< time << ",  dt "<< dt << std::endl;

    for(pHydroSolver f: schnek::BlockContainer<HydroSolver>::childBlocks())
    {
      f->timeStep(dt);
    }

    time += dt;
    ++timestep;
  }

  schnek::DiagnosticManager::instance().execute();
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  static const int eulerModelDim = EulerKnp<DIMENSION>::dim;
  static const int adiabaticModelDim = AdiabaticKnp<DIMENSION>::dim;

  try
  {
    schnek::BlockClasses blocks;

    blocks.registerBlock("vellamo").setClass<Vellamo>();
    blocks("Fields").setClass<HydroFields>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();
    blocks("CompressibleEuler").setClass<EulerKnp<DIMENSION>>();
    blocks("Adiabatic").setClass<AdiabaticKnp<DIMENSION>>();

    blocks("ZeroNeumannBoundaryEuler").setClass<ZeroNeumannBoundaryBlock<Field, eulerModelDim> >();
    blocks("ZeroNeumannBoundaryAdiabatic").setClass<ZeroNeumannBoundaryBlock<Field, adiabaticModelDim> >();
//    blocks("WallBoundary").setClass<WallBoundaryBlock>();

    blocks("vellamo").addChildren("Fields")("FieldDiag")("CompressibleEuler")("Adiabatic");
    blocks("CompressibleEuler").addChildren("ZeroNeumannBoundaryEuler")("WallBoundary");
    blocks("Adiabatic").addChildren("ZeroNeumannBoundaryAdiabatic");


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
