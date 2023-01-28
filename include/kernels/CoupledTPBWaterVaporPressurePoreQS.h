#ifndef COUPLEDTPBWATERVAPORPRESSUREPOREQS_H
#define COUPLEDTPBWATERVAPORPRESSUREPOREQS_H

#include "Kernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class CoupledTPBWaterVaporPressurePoreQS;
class Function;

template<>
InputParameters validParams<CoupledTPBWaterVaporPressurePoreQS>();

/**
 * Coupled Butler-Volmer kernel in sinh (empirical) form for WaterVapor concentration in pore.
 */

class CoupledTPBWaterVaporPressurePoreQS : public Kernel
{
public:
  CoupledTPBWaterVaporPressurePoreQS(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  // const Real _s0;       // exchange (equilibrium) current rate
  const Real _z;
  const Real _F;
  const Real _R;
  const Real _T;
  const Real _Kf;
  const Real _p_total;
  // const Function & _func_phi_Ni;
  // const Real _pO2_CE;
  const MaterialProperty<Real> & _s0;
  const MaterialProperty<Real> & _pO2_CE;

  unsigned int _num_phi_YSZ;
  const VariableValue & _phi_YSZ;
  unsigned int _num_phi_Ni;
  const VariableValue & _phi_Ni;
};

#endif //COUPLEDTPBWaterVaporPressurePoreQS_H
