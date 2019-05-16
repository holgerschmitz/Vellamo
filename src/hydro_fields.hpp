/*
 * em_fields.hpp
 *
 *  Created on: 20 Apr 2018
 *      Author: Holger Schmitz
 */


#ifndef VELLAMO_HYDRO_FIELDS_H
#define VELLAMO_HYDRO_FIELDS_H

#include "types.hpp"

#include <schnek/variables/blockcontainer.hpp>

class Vellamo;

/** @brief A container block for hydrodynamic fields
 *
 *  Multiple sets of fields can be defined
 */
class HydroFields : public schnek::ChildBlock<HydroFields>
{
  private:
    friend class Vellamo;
    /// The fluid mass density
    pField Rho;
    /// The fluid momentum
    schnek::Array<pField, DIMENSION> M;
    /// The internal energy of the fluid
    pField E;

    /// The setup parameter holding the reference to #Rho
    schnek::pParameter Rho_parameter;
    /// The setup parameters holding the reference to #M
    schnek::Array<schnek::pParameter, DIMENSION> M_parameters;
    /// The setup parameters holding the reference to #E
    schnek::pParameter E_parameter;

    /// The initialiser value for #Rho
    double initRho;
    /// The initialiser value for #M
    Vector initM;
    /// The initialiser value for #E
    double initE;

    void fillValues();
  public:

    /** @brief Constructor taking an optional parent block
     *
     * Constructor is compatible with the schnek::Block constructor
     *
     */
    HydroFields(schnek::pBlock parent = schnek::pBlock()) : schnek::ChildBlock<HydroFields>(parent)
    {}

    /** @brief Register the hydrodynamic fields
     *
     * The hydrodynamic fields are created and registered with the Block storage.
     */
    void registerData();
    void initParameters(schnek::BlockParameters &parameters);
    void init();
};

#endif
