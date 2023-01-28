#ifndef COUPLEDTPBPOTENTIALYSZQS_H
#define COUPLEDTPBPOTENTIALYSZQS_H

#include "Kernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class CoupledTPBPotentialYSZQS;
class Function;

template<>
InputParameters validParams<CoupledTPBPotentialYSZQS>();

/**
 * Coupled Butler-Volmer TPB kernel in sinh() form for potential in YSZ.
 */

class CoupledTPBPotentialYSZQS : public Kernel
{
public:
  CoupledTPBPotentialYSZQS(const InputParameters & parameters);

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

  unsigned int _num_p_O2;
  const VariableValue & _p_O2;
  unsigned int _num_phi_LSM;
  const VariableValue & _phi_LSM;
};

#endif //COUPLEDTPBPOTENTIALYSZQS_H
