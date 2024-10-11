//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IsotropicPlasticityTempDepStressUpdate.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", ADIsotropicPlasticityTempDepStressUpdate);
registerMooseObject("SolidMechanicsApp", IsotropicPlasticityTempDepStressUpdate);

template <bool is_ad>
InputParameters
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the discrete material in a radial return isotropic "
                             "plasticity model.  This class is one of the basic radial return "
                             "constitutive models, yet it can be used in conjunction with other "
                             "creep and plasticity materials for more complex simulations.");
  // Linear strain hardening parameters
  params.addParam<FunctionName>("yield_stress_function",
                                "Yield stress as a function of temperature");
  params.addParam<Real>("yield_stress", "The point at which plastic strain begins accumulating");
  // params.addParam<FunctionName>("hardening_function",
  //                               "True stress as a function of plastic strain");
  params.addParam<std::vector<FunctionName>>("hardening_function_list",
                                "Vectors of stress vs strain hardening slopes");
  params.addParam<std::vector<Real>>("hardening_temps",
                                "Temperatures for hardening slopes");
  params.addParam<Real>("hardening_constant", "Hardening slope");
  params.addCoupledVar("temperature", 0.0, "Coupled Temperature");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  return params;
}

template <bool is_ad>
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::IsotropicPlasticityTempDepStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),
    _plastic_prepend(this->template getParam<std::string>("plastic_prepend")),
    _yield_stress_function(this->isParamValid("yield_stress_function")
                               ? &this->getFunction("yield_stress_function")
                               : nullptr),
    _yield_stress(this->isParamValid("yield_stress") ? this->template getParam<Real>("yield_stress")
                                                     : 0),
    _hardening_constant(this->isParamValid("hardening_constant")
                            ? this->template getParam<Real>("hardening_constant")
                            : 0),
    // _hardening_function(this->isParamValid("hardening_function")
    //                         ? &this->getFunction("hardening_function")
    //                         : nullptr), 
    _hardening_function_list(this->template getParam<std::vector<FunctionName>>("hardening_function_list")),
    _hardening_temps(this->template getParam<std::vector<Real>>("hardening_temps")),
    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _hardening_slope(0.0),
    _plastic_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>(
        _base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>(
        _base_name + _plastic_prepend + "plastic_strain")),
    _hardening_variable(
        this->template declareGenericProperty<Real, is_ad>(_base_name + "hardening_variable")),
    _hardening_variable_old(
        this->template getMaterialPropertyOld<Real>(_base_name + "hardening_variable")),
    _temperature(this->template coupledGenericValue<is_ad>("temperature"))
{
  if (parameters.isParamSetByUser("yield_stress") && _yield_stress <= 0.0)
    mooseError("Yield stress must be greater than zero");

  // Both of these parameters are given default values by derived classes, which makes them valid
  if (_yield_stress_function == nullptr && !this->isParamValid("yield_stress"))
    mooseError("Either yield_stress or yield_stress_function must be given");
  // if (!parameters.isParamValid("hardening_constant") && !this->isParamValid("hardening_function"))
  //   mooseError("Either hardening_constant or hardening_function must be defined");

  // if (parameters.isParamSetByUser("hardening_constant") && this->isParamValid("hardening_function"))
  //   mooseError(
  //       "Only the hardening_constant or only the hardening_function can be defined but not both");
  unsigned int num_hardening_funcs = _hardening_function_list.size();
  unsigned int num_hardening_temps = _hardening_temps.size();
  if (num_hardening_funcs != num_hardening_temps)
    mooseError("Number of hardening_function_list (",
               num_hardening_funcs,
               ") must match the number of hardening_temps (",
               num_hardening_temps,
               ") for the hardening slopes.");
  
  _hardening_functions.resize(num_hardening_funcs);
  for (unsigned int i = 0; i < num_hardening_funcs; i++)
  {
    _hardening_functions[i] = &this->getFunctionByName(_hardening_function_list[i]);
  }
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::initQpStatefulProperties()
{
  _hardening_variable[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  computeYieldStress(elasticity_tensor);

  _yield_condition = effective_trial_stress - _hardening_variable_old[_qp] - _yield_stress;
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeResidual(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {
    _hardening_slope = computeHardeningDerivative(scalar);
    _hardening_variable[_qp] = computeHardeningValue(scalar);

    return (effective_trial_stress - _hardening_variable[_qp] - _yield_stress) /
               _three_shear_modulus -
           scalar;
  }

  return 0.0;
}

template <bool is_ad>
GenericReal<is_ad>
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /*effective_trial_stress*/, const GenericReal<is_ad> & /*scalar*/)
{
  if (_yield_condition > 0.0)
    return -1.0 - _hardening_slope / _three_shear_modulus;

  return 1.0;
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::iterationFinalize(const GenericReal<is_ad> & scalar)
{
  if (_yield_condition > 0.0)
    _hardening_variable[_qp] = computeHardeningValue(scalar);
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
}

template <bool is_ad>
GenericReal<is_ad>
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeHardeningValue(
    const GenericReal<is_ad> & scalar)
{
  // if (_hardening_function)
  // {
  //   const Real strain_old = this->_effective_inelastic_strain_old[_qp];
  //   return _hardening_function->value(strain_old + scalar) - _yield_stress;
  // }
  
  if (_hardening_functions[0])
  {
    const Real strain_old = this->_effective_inelastic_strain_old[_qp];
    if (_temperature[_qp] < _hardening_temps[0])
      return _hardening_functions[0]->value(strain_old + scalar) - _yield_stress;
    else if (_temperature[_qp] >= _hardening_temps[-1])
      return _hardening_functions[-1]->value(strain_old + scalar) - _yield_stress;
  }


  return _hardening_variable_old[_qp] + _hardening_slope * scalar;
}

template <bool is_ad>
GenericReal<is_ad>
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeHardeningDerivative(
    const GenericReal<is_ad> & /*scalar*/)
{
  // if (_hardening_function)
  // {
  //   const Real strain_old = this->_effective_inelastic_strain_old[_qp];
  //   return _hardening_function->timeDerivative(strain_old);
  // }
  if (_hardening_functions[0])
  {
    const Real strain_old = this->_effective_inelastic_strain_old[_qp];
    if (_temperature[_qp] < _hardening_temps[0])
      return _hardening_functions[0]->timeDerivative(strain_old);
    else if (_temperature[_qp] >= _hardening_temps[-1])
      return _hardening_functions[-1]->timeDerivative(strain_old);
  }

  return _hardening_constant;
}

template <bool is_ad>
void
IsotropicPlasticityTempDepStressUpdateTempl<is_ad>::computeYieldStress(
    const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/)
{
  if (_yield_stress_function)
  {
    static const Moose::GenericType<Point, is_ad> p;
    _yield_stress = _yield_stress_function->value(_temperature[_qp], p);

    if (_yield_stress <= 0.0)
      mooseError("In ",
                 this->_name,
                 ": The calculated yield stress (",
                 _yield_stress,
                 ") is less than zero");
  }
}

template class IsotropicPlasticityTempDepStressUpdateTempl<false>;
template class IsotropicPlasticityTempDepStressUpdateTempl<true>;
