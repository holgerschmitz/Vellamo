/*
 * mpulse.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#include "diagnostic.hpp"

#include <schnek/parser.hpp>
#include <schnek/tools/fieldtools.hpp>
#include <schnek/tools/literature.hpp>
#include <schnek/util/logger.hpp>

#include <boost/foreach.hpp>

#include <mpi.h>

#include <fstream>
#include <string>
#include <unistd.h>

Vellamo::Vellamo()
{
  instance = this;
}

void Vellamo::initParameters(schnek::BlockParameters &parameters)
{
  parameters.addArrayParameter("N", gridSize, 100);
  parameters.addArrayParameter("L", size);
  parameters.addParameter("tMax", &tMax);
  parameters.addParameter("cflFactor", &cflFactor, 0.99);
  x_parameters = parameters.addArrayParameter("", x, schnek::BlockParameters::readonly);

  Rho_parameter = parameters.addParameter("rho", &initRho, 0.0);

  V_parameters = parameters.addArrayParameter("Vx", initV, 0.0);

  E_parameter = parameters.addParameter("E", &initE, 0.0);

  spaceVars = schnek::pParametersGroup(new schnek::ParametersGroup());
  spaceVars->addArray(x_parameters);
}

void Vellamo::registerData()
{
  addData("Rho", Rho);

  addData("Vx", Vx);
  addData("Vy", Vy);
  addData("Vz", Vz);

  addData("E", E);
}


void Vellamo::fillValues()
{
  schnek::pBlockVariables blockVars = getVariables();
  schnek::pDependencyMap depMap(new schnek::DependencyMap(blockVars));

  schnek::DependencyUpdater updater(depMap);

  updater.addIndependentArray(x_parameters);
  schnek::fill_field(*Rho, x, initRho, updater, Rho_parameter);

  schnek::fill_field(*Vx, x, initV[0], updater, V_parameters[0]);
  schnek::fill_field(*Vy, x, initV[1], updater, V_parameters[1]);
  schnek::fill_field(*Vz, x, initV[2], updater, V_parameters[2]);

  schnek::fill_field(*E, x, initE, updater, E_parameter);
}

void Vellamo::init()
{
  globalMax = gridSize - 2;

  subdivision.init(gridSize, 2);

  dx = size / gridSize;
  dt = cflFactor*std::min(dx[0],std::min(dx[1],dx[2]))/clight;

  Index low  = subdivision.getLo();
  Index high = subdivision.getHi();

  Index lowIn  = subdivision.getInnerLo();
  Index highIn = subdivision.getInnerHi();

  innerRange = Range(lowIn, highIn);
  schnek::Range<double, DIMENSION> domainSize(schnek::Array<double, DIMENSION>(lowIn[0],lowIn[1]), schnek::Array<double, DIMENSION>(highIn[0],highIn[1]));
  schnek::Array<bool, DIMENSION> stagger;

  stagger = false;

  Rho = new Field(lowIn, highIn, domainSize, stagger, 2);

  Vx = new Field(lowIn, highIn, domainSize, stagger, 2);
  Vy = new Field(lowIn, highIn, domainSize, stagger, 2);
  Vz = new Field(lowIn, highIn, domainSize, stagger, 2);

  E = new Field(lowIn, highIn, domainSize, stagger, 2);

  fillValues();
}

void Vellamo::execute()
{
  int step = 0;
  double time = 0.0;


  while (time<=tMax)
  {
    schnek::DiagnosticManager::instance().execute();

    if (subdivision.master())
      schnek::Logger::instance().out() <<"Time "<< time << std::endl;

      BOOST_FOREACH(Solver *f, getChildren())
      {
        f->rungeKuttaStep(dt);
      }
  }

  schnek::DiagnosticManager::instance().execute();
}

int main (int argc, char** argv) {

  MPI_Init(&argc, &argv);

  try
  {
    schnek::BlockClasses blocks;

    blocks.registerBlock("mpulse").setClass<Vellamo>();
    blocks("FieldDiag").setClass<FieldDiagnostic>();

    blocks("mpulse").addChildren("FieldDiag");

    std::ifstream in("mpulse.setup");
    if (!in) throw std::string("Could not open file 'mpulse.setup'");

    schnek::Parser P("mpulse", "mpulse", blocks);
    registerCMath(P.getFunctionRegistry());

    //P.getFunctionRegistry().registerFunction("random",randomRange);

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
