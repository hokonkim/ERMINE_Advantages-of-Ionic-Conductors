//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "MaterialProperty.h"

class VolCr_prev : public Material
{
public:
  VolCr_prev(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  // virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

private:
  Real _initial_VolCr_prev;
  const MaterialProperty<Real> & _VolCr;
  /**
   * Create two MooseArray Refs to hold the current
   * and previous material properties respectively
   */
  MaterialProperty<Real> & _VolCr_prev;
  // const MaterialProperty<Real> & _VolCr_old;
};
