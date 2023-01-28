#ifndef COUPLEDTPBOXYGENPRESSUREPOREQS_H
#define COUPLEDTPBOXYGENPRESSUREPOREQS_H

#include "Kernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class CoupledTPBOxygenPressurePoreQS;
class Function;

template<>
InputParameters validParams<CoupledTPBOxygenPressurePoreQS>();

/**
 * Coupled Butler-Volmer kernel in sinh (empirical) form for oxygen concentration in pore.
 */

class CoupledTPBOxygenPressurePoreQS : public Kernel
{
public:
  CoupledTPBOxygenPressurePoreQS(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  // const Real _s0;
  const Real _z;
  const Real _F;
  const Real _R;
  const Real _T;
  // const Function & _func_phi_LSM;
  // const Real _pO2_CE;
  const MaterialProperty<Real> & _s0;  // exchange (equilibrium) current rate
  const MaterialProperty<Real> & _pO2_CE;

  unsigned int _num_phi_YSZ;
  const VariableValue & _phi_YSZ;
  unsigned int _num_phi_LSM;
  const VariableValue & _phi_LSM;
};

#endif //COUPLEDTPBOXYGENPRESSUREPOREQS_H
