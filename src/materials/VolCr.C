//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolCr.h"

registerMooseObject("ermineApp", VolCr);

InputParameters
VolCr::validParams()
{
  InputParameters params = Material::validParams();

  // Allow users to specify vectors defining the points of a piecewise function formed via linear
  // interpolation.
  params.addParam<Real>("initial_VolCr", "The initial volume of Cr at tpb");
  params.addParam<MaterialPropertyName>("VolCr_new", "The amount of Cr to be added");
  return params;
}

VolCr::VolCr(const InputParameters & parameters)
  : Material(parameters),
  _initial_VolCr(getParam<Real>("initial_VolCr")),
  // Declare that this material is going to have a Real
  // valued property named "VolCr" that Kernels can use.
  _VolCr_new(getMaterialProperty<Real>("VolCr_new")),
  // Get a parameter value for the VolCr
  _VolCr(declareProperty<Real>("VolCr")),
  // Retrieve/use an old value of VolCr.
  // Note: this is _expensive_ as we have to store values for each
  // qp throughout the mesh. Only do this if you REALLY need it!
  _VolCr_old(getMaterialPropertyOld<Real>("VolCr"))
{
}

void
VolCr::initQpStatefulProperties()
{
  // init the VolCr property (this will become
  // _VolCr_old in the first call of computeProperties)
  _VolCr[_qp] = _initial_VolCr;
}

void
VolCr::computeQpProperties()
{
  // VolCr is obtained from VolCr_old and VolCr_new
  _VolCr[_qp] = _VolCr_old[_qp] + _VolCr_new[_qp];
}
