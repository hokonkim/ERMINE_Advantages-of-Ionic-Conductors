//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolCr_prev.h"

registerMooseObject("ermineApp", VolCr_prev);

InputParameters
VolCr_prev::validParams()
{
  InputParameters params = Material::validParams();

  // Allow users to specify vectors defining the points of a piecewise function formed via linear
  // interpolation.
  params.addParam<Real>("initial_VolCr_prev", "The initial volume of Cr at tpb");
  params.addParam<MaterialPropertyName>("VolCr", "The amount of Cr previously at tpb");
  return params;
}

VolCr_prev::VolCr_prev(const InputParameters & parameters)
  : Material(parameters),
  _initial_VolCr_prev(getParam<Real>("initial_VolCr_prev")),
  // Declare that this material is going to have a Real
  // valued property named "VolCr" that Kernels can use.
  _VolCr(getMaterialPropertyOld<Real>("VolCr")),
  // Get a parameter value for the VolCr
  _VolCr_prev(declareProperty<Real>("VolCr_prev"))
  // Retrieve/use an old value of VolCr.
  // Note: this is _expensive_ as we have to store values for each
  // qp throughout the mesh. Only do this if you REALLY need it!
  // _VolCr_old(getMaterialPropertyOld<Real>("VolCr_prev"))
{
}

// void
// VolCr_prev::initQpStatefulProperties()
// {
//   // init the VolCr property (this will become
//   // _VolCr_old in the first call of computeProperties)
//   _VolCr_prev[_qp] = _initial_VolCr_prev;
// }

void
VolCr_prev::computeQpProperties()
{
  // VolCr_prev is obtained from VolCr
  _VolCr_prev[_qp] = _VolCr[_qp];
}
