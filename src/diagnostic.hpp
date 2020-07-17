/*
 * diagnostic.hpp
 *
 *  Created on: 5 Oct 2016
 *      Author: Holger Schmitz
 */

#ifndef DIAGNOSTIC_HPP_
#define DIAGNOSTIC_HPP_

#include "vellamo.hpp"
#include <schnek/diagnostic/diagnostic.hpp>
#include <schnek/diagnostic/hdfdiagnostic.hpp>

class FieldDiagnostic : public schnek::HDFGridDiagnostic<Field, pField, schnek::DeltaTimeDiagnostic>
{
  protected:
    Index getGlobalMin() override;
    Index getGlobalMax() override;
//    schnek::pHdfAttributes getAttributes();
};



#endif /* DIAGNOSTIC_HPP_ */
