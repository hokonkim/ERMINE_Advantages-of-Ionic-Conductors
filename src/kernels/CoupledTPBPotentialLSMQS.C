#include "CoupledTPBPotentialLSMQS.h"
#include "Function.h"
#include <cmath>

registerMooseObject("ermineApp", CoupledTPBPotentialLSMQS);


template<>
InputParameters validParams<CoupledTPBPotentialLSMQS>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Fully Coupled TPB reaction source kernel in Butler-Volmer sinh() form for potential in LSM. j = 2 * j0 * sinh(coef * eta)");
  // params.addRequiredParam<Real>("s0", "Exchange volumetric current density rate (A/cm^3)");
  params.addParam<Real>("z", 4.0, "electron number (num of electrons transferred)");
  params.addParam<Real>("F", 96485.33289, "Faraday constant (C/mol)");
  params.addParam<Real>("R", 8.3144598, "Gas constant (J/K/mol)");
  params.addParam<Real>("T", 1073.0, "Temperature (T)");
  // params.addRequiredParam<FunctionName>("function_phi_LSM", "Function for the LSM potential");
  // params.addRequiredParam<Real>("pO2_CE", "Oxygen partial pressure at the counter electrode (atm)");
  params.addRequiredParam<MaterialPropertyName>("s0", "Exchange volumetric current density rate (A/cm^3)");
  params.addRequiredParam<MaterialPropertyName>("pO2_CE", "Oxygen partial pressure at the counter electrode (atm) (material property)");
  params.addRequiredCoupledVar("p_O2", "The coupled oxygen partial pressure variable in pore");
  params.addRequiredCoupledVar("phi_YSZ", "The coupled potential variable in eta = potential_rev - (phi_LSM - phi_YSZ)");
  return params;
}

CoupledTPBPotentialLSMQS::CoupledTPBPotentialLSMQS(const InputParameters & parameters) :
    Kernel(parameters),
    // _s0(getParam<Real>("s0")),
    _z(getParam<Real>("z")),
    _F(getParam<Real>("F")),
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    // _func_phi_LSM(getFunction("function_phi_LSM")),
    // _pO2_CE(getParam<Real>("pO2_CE")),
    _s0(getMaterialProperty<Real>("s0")),
    _pO2_CE(getMaterialProperty<Real>("pO2_CE")),
    _num_p_O2(coupled("p_O2")),
    _p_O2(coupledValue("p_O2")),
    _num_phi_YSZ(coupled("phi_YSZ")),
    _phi_YSZ(coupledValue("phi_YSZ"))
{
}

// Weak form is (somthing) = 0
// What I need to implement is the term in (something)

Real
CoupledTPBPotentialLSMQS::computeQpResidual()
{
  Real b       = _z * _F / _R / _T;
  Real E_conc  = - _R * _T / 4 / _F * log(_pO2_CE[_qp] / _p_O2[_qp]);
  // Real phi_LSM = _func_phi_LSM.value(_t, _q_point[_qp]);
  Real eta     = E_conc - (_u[_qp] - _phi_YSZ[_qp]);

  Real res     = -2 * _s0[_qp] * sinh(0.5 * b * eta);
  return res * _test[_i][_qp];
}

// Jacobian is the derivative of residual with respect to other coupled variables

Real
CoupledTPBPotentialLSMQS::computeQpJacobian()
{
  Real b         = _z * _F / _R / _T;
  Real E_conc    = - _R * _T / 4 / _F * log(_pO2_CE[_qp] / _p_O2[_qp]);
  // Real phi_LSM   = _func_phi_LSM.value(_t, _q_point[_qp]);
  Real eta     = E_conc - (_u[_qp] - _phi_YSZ[_qp]);

  Real ds_deta   = b * _s0[_qp] * cosh(0.5 * b * eta);
  Real deta_dLSM = -1;

  Real jac       = -ds_deta * deta_dLSM;
  return jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
CoupledTPBPotentialLSMQS::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _num_p_O2)
  {
    Real b       = _z * _F / _R / _T;
    Real E_conc  = - _R * _T / 4 / _F * log(_pO2_CE[_qp] / _p_O2[_qp]);
    // Real phi_LSM = _func_phi_LSM.value(_t, _q_point[_qp]);
    Real eta     = E_conc - (_u[_qp] - _phi_YSZ[_qp]);

    Real ds_deta = b * _s0[_qp] * cosh(0.5 * b * eta);
    Real deta_dE = 1;
    Real dE_dpO2 = _R * _T / 4 / _F / _p_O2[_qp];

    Real jac     = -ds_deta * deta_dE * dE_dpO2;
    return jac * _test[_i][_qp] * _phi[_j][_qp];
  }
  else if (jvar == _num_phi_YSZ)
  {
    Real b       = _z * _F / _R / _T;
    Real E_conc  = - _R * _T / 4 / _F * log(_pO2_CE[_qp] / _p_O2[_qp]);
    Real eta     = E_conc - (_u[_qp] - _phi_YSZ[_qp]);

    Real ds_deta   = b * _s0[_qp] * cosh(0.5 * b * eta);
    Real deta_dphi = 1;

    Real jac     = -ds_deta * deta_dphi;
    return jac * _test[_i][_qp] * _phi[_j][_qp];
  }
  else return 0.0;
}
