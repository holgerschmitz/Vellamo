/*
 * mpulse.hpp
 *
 *  Created on: 4 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef MPULSE_MPULSE_H
#define MPULSE_MPULSE_H

#include <schnek/grid.hpp>
#include <schnek/variables.hpp>

#ifdef NDEBUG
#define MPulseGridChecker schnek::GridNoArgCheck
#else
#define MPulseGridChecker schnek::GridAssertCheck
#endif

static const size_t DIMENSION = 2;

typedef schnek::Array<int, DIMENSION> Index;
typedef schnek::Array<double, DIMENSION> Vector;
typedef schnek::Field<double, DIMENSION> Field;
typedef boost::shared_ptr<Field> pField;
typedef schnek::Range<int, DIMENSION> Range;
typedef schnek::Array<bool, DIMENSION> Stagger;

static const double clight = 299792458;
static const double clight2 = clight*clight;

class Solver;

class Vellamo : public schnek::Block, schnek::BlockContainer<Solver>
{
  private:
    static Vellamo *instance;
    Index globalMax;
    Index gridSize;
    Vector size;
    Vector dx;

    double cflFactor;
    double dt;


    double tMax;
    Range innerRange;
    pField Rho;
    schnek::Array<pField, DIMENSION> M;
    pField E;

    schnek::MPICartSubdivision<Field> subdivision;

    Vector x;
    schnek::Array<schnek::pParameter, DIMENSION> x_parameters;
    schnek::pParameter Rho_parameter;
    schnek::Array<schnek::pParameter, DIMENSION> M_parameters;
    schnek::pParameter E_parameter;
    schnek::pParametersGroup spaceVars;

    double initRho;
    Vector initV;
    double initE;
  protected:
    void initParameters(schnek::BlockParameters &blockPars);
    void registerData();
    void fillValues();
  public:
    Vellamo();
    void init();
    void execute();

    static Index getGlobalMax() { return instance->globalMax; }
    static Vector getDx() { return instance->dx; }
    static schnek::DomainSubdivision<Field> &getSubdivision() { return instance->subdivision; };
};

#endif // MPULSE_MPULSE_H
