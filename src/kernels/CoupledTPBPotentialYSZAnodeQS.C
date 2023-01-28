#include "CoupledTPBPotentialYSZAnodeQS.h"
#include "Function.h"
#include <cmath>

registerMooseObject("ermineApp", CoupledTPBPotentialYSZAnodeQS);


template<>
InputParameters validParams<CoupledTPBPotentialYSZAnodeQS>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Fully Coupled TPB reaction source kernel in Butler-Volmer sinh() form for potential in YSZ. j = 2 * j0 * sinh(coef * eta)");
  // params.addRequiredParam<Real>("s0", "Exchange volumetric current density rate (A/cm^3)");
  params.addParam<Real>("z", 4.0, "Electron number (num of electrons transferred)");
  params.addParam<Real>("F", 96485.33289, "Faraday constant (C/mol)");
  params.addParam<Real>("R", 8.3144598, "Gas constant (J/K/mol)");
  params.addParam<Real>("T", 1073.0, "Temperature (T)");
  params.addParam<Real>("Kf", 7.638e8, "Equilibrium constant for H2 + 1/2O2 = H2O at 1073 K");
  params.addRequiredParam<Real>("p_total", "Total pressure of Hydrogen and Water Vapor at the anode (atm)");
  // params.addRequiredParam<FunctionName>("function_phi_Ni", "Function for the Ni potential");
  // params.addRequiredParam<Real>("pO2_CE", "Oxygen partial pressure at the counter electrode (atm)");
  params.addRequiredParam<MaterialPropertyName>("s0", "Exchange volumetric current density rate (A/cm^3)");
  params.addRequiredParam<MaterialPropertyName>("pO2_CE", "Oxygen partial pressure at the counter electrode (atm) (material property)");
  params.addRequiredCoupledVar("p_H2O", "The coupled Water Vapor partial pressure variable in pore");
  params.addRequiredCoupledVar("phi_Ni", "The coupled potential variable in eta = potential_rev - (phi_YSZ - phi_Ni)");
  return params;
}

CoupledTPBPotentialYSZAnodeQS::CoupledTPBPotentialYSZAnodeQS(const InputParameters & parameters) :
    Kernel(parameters),
    // _s0(getParam<Real>("s0")),
    _z(getParam<Real>("z")),
    _F(getParam<Real>("F")),
    _R(getParam<Real>("R")),
    _T(getParam<Real>("T")),
    _Kf(getParam<Real>("Kf")),
    _p_total(getParam<Real>("p_total")),
    // _func_phi_Ni(getFunction("function_phi_Ni")),
    // _pO2_CE(getParam<Real>("pO2_CE")),
    _s0(getMaterialProperty<Real>("s0")),
    _pO2_CE(getMaterialProperty<Real>("pO2_CE")),
    _num_p_H2O(coupled("p_H2O")),
    _p_H2O(coupledValue("p_H2O")),
    _num_phi_Ni(coupled("phi_Ni")),
    _phi_Ni(coupledValue("phi_Ni"))
{
}

// Weak form is (somthing) = 0
// What I need to implement is the term in (something)

Real
CoupledTPBPotentialYSZAnodeQS::computeQpResidual()
{
  Real b       = _z * _F / _R / _T;
  Real E_conc  = _R * _T / 4 / _F * log(_pO2_CE[_qp] / (pow (_p_H2O[_qp]/((_p_total-_p_H2O[_qp]) * _Kf), 2.0)));
  // Real phi_Ni  = _func_phi_Ni.value(_t, _q_point[_qp]);
  Real eta     = E_conc - (_u[_qp] - _phi_Ni[_qp]);

  Real res     = -2 * _s0[_qp] * sinh(0.5 * b * eta); // it doesn't need to be changed, only z should be changed.
  return res * _test[_i][_qp];
}

// Jacobian is the derivative of residual with respect to other coupled variables

Real
CoupledTPBPotentialYSZAnodeQS::computeQpJacobian()
{
  Real b       = _z * _F / _R / _T;
  Real E_conc  = _R * _T / 4 / _F * log(_pO2_CE[_qp] / (pow (_p_H2O[_qp]/((_p_total-_p_H2O[_qp]) * _Kf), 2.0)));
  // Real phi_Ni  = _func_phi_Ni.value(_t, _q_point[_qp]);
  Real eta     = E_conc - (_u[_qp] - _phi_Ni[_qp]);

  Real ds_deta   = b * _s0[_qp] * cosh(0.5 * b * eta);
  Real deta_dphi = -1;

  Real jac       = -ds_deta * deta_dphi; // it doesn't need to be changed, only z should be changed.
  return jac * _test[_i][_qp] * _phi[_j][_qp];
}

Real
CoupledTPBPotentialYSZAnodeQS::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _num_p_H2O)
  {
    Real b       = _z * _F / _R / _T;
    Real E_conc  = _R * _T / 4 / _F * log(_pO2_CE[_qp] / (pow (_p_H2O[_qp]/((_p_total-_p_H2O[_qp]) * _Kf), 2.0)));
    // Real phi_Ni  = _func_phi_Ni.value(_t, _q_point[_qp]);
    Real eta     = E_conc - (_u[_qp] - _phi_Ni[_qp]);

    Real ds_deta = b * _s0[_qp] * cosh(0.5 * b * eta);
    Real deta_dE = 1;
    Real dE_dpO2 = -1 / (b * (pow (_p_H2O[_qp]/((_p_total - _p_H2O[_qp]) * _Kf), 2.0)));
    Real dpO2_dpH2O = ((2 * _p_H2O[_qp]) / (pow(_p_total - _p_H2O[_qp], 2.0) * pow(_Kf, 2.0))) * (1 + (_p_H2O[_qp] / (_p_total - _p_H2O[_qp])));

    Real jac     = -ds_deta * deta_dE * dE_dpO2 * dpO2_dpH2O; // it doesn't need to be changed, only z should be changed.
    return jac * _test[_i][_qp] * _phi[_j][_qp];
  }
  else if (jvar == _num_phi_Ni)
  {
    Real b       = _z * _F / _R / _T;
    Real E_conc  = _R * _T / 4 / _F * log(_pO2_CE[_qp] / (pow (_p_H2O[_qp]/((_p_total-_p_H2O[_qp]) * _Kf), 2.0)));
    // Real phi_Ni  = _func_phi_Ni.value(_t, _q_point[_qp]);
    Real eta     = E_conc - (_u[_qp] - _phi_Ni[_qp]);

    Real ds_deta  = b * _s0[_qp] * cosh(0.5 * b * eta);
    Real deta_dNi = 1;

    Real jac     = -ds_deta * deta_dNi;
    return jac * _test[_i][_qp] * _phi[_j][_qp];
  }
  else return 0.0;
}
