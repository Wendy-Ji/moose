//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeMeanThermalExpansionFunctionPhaseEigenstrain.h"
#include "Function.h"
#include "PiecewiseLinear.h"

registerMooseObject("SolidMechanicsApp", ComputeMeanThermalExpansionFunctionPhaseEigenstrain);
registerMooseObject("SolidMechanicsApp", ADComputeMeanThermalExpansionFunctionPhaseEigenstrain);

template <bool is_ad>
InputParameters
ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<is_ad>::validParams()
{
  InputParameters params = ComputeMeanThermalExpansionEigenstrainBaseTempl<is_ad>::validParams();
  params.addClassDescription("Computes eigenstrain due to thermal expansion using a function that "
                             "describes the mean thermal expansion as a function of temperature");
  params.addRequiredParam<std::vector<FunctionName>>(
      "thermal_expansion_functions",
      "Vector of functions describing the mean thermal expansion as a function of temperature, in order of phase (aust, bain, mart, parent, ferrite)");
  params.addRequiredParam<Real>("thermal_expansion_function_reference_temperature",
                                "Reference temperature for thermal_exansion_function (IMPORTANT: "
                                "this is different in general from the stress_free_temperature)");
  // params.addCoupledVar("phase", "Array variable for phase proportions in order of aust, bain, mart, parent, ferrite");
  return params;
}

template <bool is_ad>
ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<
    is_ad>::ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl(const InputParameters & parameters)
  : ComputeMeanThermalExpansionEigenstrainBaseTempl<is_ad>(parameters),
    _thermal_expansion_functions_names(
        this->template getParam<std::vector<FunctionName>>("thermal_expansion_functions")),
    _thexp_func_ref_temp(
        this->template getParam<Real>("thermal_expansion_function_reference_temperature"))
    // _phase(this->template getParam<coupledArrayValue>("phase"))
{
  const unsigned int len = _thermal_expansion_functions_names.size();
  for (unsigned int i = 1; i < len; ++i)
  {
    const PiecewiseLinear * const f = dynamic_cast<const PiecewiseLinear *>(
    &this->getFunctionByName(_thermal_expansion_functions_names[i]));
    _thermal_expansion_functions[i] = f;
  }
}

template <bool is_ad>
Real
ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<is_ad>::referenceTemperature()
{
  return _thexp_func_ref_temp;
}

template <bool is_ad>
ValueAndDerivative<is_ad>
ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<is_ad>::meanThermalExpansionCoefficient(
    const ValueAndDerivative<is_ad> & temperature)
{
  // // we need these two branches because we cannot yet evaluate Functions with ChainedReals
  // if constexpr (is_ad)
  //   return _thermal_expansion_function.value(temperature);
  // else
  //   return {_thermal_expansion_function.value(temperature.value()),
  //           _thermal_expansion_function.timeDerivative(temperature.value()) *
  //               temperature.derivatives()};

  // return _thermal_expansion_functions[0].value(temperature)*_phase[_qp][0] +_thermal_expansion_functions[1].value(temperature)*_phase[_qp][1] + _thermal_expansion_functions[2].value(temperature)*_phase[_qp][2] + _thermal_expansion_functions[3].value(temperature)*_phase[_qp][3] + _thermal_expansion_functions[4].value(temperature)*_phase[_qp][4];
  return _thermal_expansion_functions[0]->value(temperature,Point())*0;
}

template class ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<false>;
template class ComputeMeanThermalExpansionFunctionPhaseEigenstrainTempl<true>;

