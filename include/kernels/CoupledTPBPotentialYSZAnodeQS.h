#ifndef COUPLEDTPBPOTENTIALYSZANODEQS_H
#define COUPLEDTPBPOTENTIALYSZANODEQS_H

#include "Kernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class CoupledTPBPotentialYSZAnodeQS;
class Function;

template<>
InputParameters validParams<CoupledTPBPotentialYSZAnodeQS>();

/**
 * Coupled Butler-Volmer TPB kernel in sinh() form for potential in YSZ.
 */

class CoupledTPBPotentialYSZAnodeQS : public Kernel
{
public:
  CoupledTPBPotentialYSZAnodeQS(const InputParameters & parameters);

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

  unsigned int _num_p_H2O;
  const VariableValue & _p_H2O;
  unsigned int _num_phi_Ni;
  const VariableValue & _phi_Ni;
};

#endif //CoupledTPBPotentialYSZAnodeQS_H
