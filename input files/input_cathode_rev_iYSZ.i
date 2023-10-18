[Mesh]
  [./file]
    type = FileMeshGenerator
    # file = /Users/HOKON/projects/ermine/inputs/tpb_gan0002.inp
    file = /Users/HOKON/projects/ermine/inputs/test_cyl.inp
  [../]

  [./meshScale]
    type = TransformGenerator
    input = file
    transform = SCALE
    vector_value = '1e-1 1e-1 1e-1' # from mm to cm
  [../]

  [./phase_12_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = meshScale
    # master_block = 'PT_PORE_TET4'
    primary_block = 'PT_PORE_TET4'
    paired_block = 'PT_LSM_TET4'
    # paired_block = 'PT_NI_TET4'
    new_boundary = 'interface_12'
  [../]

  [./phase_13_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = phase_12_interface
    # master_block = 'PT_PORE_TET4'
    primary_block = 'PT_PORE_TET4'
    paired_block = 'PT_YSZ_TET4'
    new_boundary = 'interface_13'
  [../]

  [./phase_23_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = phase_13_interface
    # master_block = 'PT_LSM_TET4'
    primary_block = 'PT_LSM_TET4'
    # primary_block = 'PT_NI_TET4'
    paired_block = 'PT_YSZ_TET4'
    new_boundary = 'interface_23'
  [../]

  [./phase_16_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = phase_23_interface
    master_block = 'PT_PORE_TET4'
    paired_block = 'PT_IYSZ_TET4'
    new_boundary = 'interface_16'
  [../]

  [./phase_26_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = phase_16_interface
    master_block = 'PT_LSM_TET4'
    paired_block = 'PT_IYSZ_TET4'
    new_boundary = 'interface_26'
  [../]
[]

#==========================================================================#

[GlobalParams]
  pO2_CE  = 7e-21   # (atm)
  # function_phi_LSM = 'funcPotentialLSM' # (V)
  # s0 = 1e3 # (A/cm^3) (6.8 * 20)
[../]

#==========================================================================#

[Functions]
  [./funcPotentialLSM]
    type = ParsedFunction
    value = 'E_rev + E_YSZ - (eta * t)'
    vars  = 'E_rev   E_YSZ   eta'
    vals  = '1.03665 0.03607 0.4'
  [../]

  [./funcOverpotential]
    type = ParsedFunction
    value = 'eta * t'
    vars  = 'eta'
    vals  = '0.4'
  [../]

  [./funcTimeStepper]
    type = PiecewiseLinear
    x = '0.0  0.1 0.5 1.0' # time_t
    y = '0.01 0.1 0.1 0.2' # time_dt
  [../]
[]

#==========================================================================#

[Materials]
  # PT_PORE_TET4 = pore
  # PT_LSM_TET4 = LSM
  # PT_YSZ_TET4 = YSZ
  # PT_TPB_TET4 = TPB
  [./gasDiffFluxCoefficient]
    type  = ParsedMaterial
    block = 'PT_PORE_TET4 PT_TPB_TET4'
    f_name = 'gasDiffFluxCoef'
    constant_names        = 'D_O2  R          T'        # (cm^2/s), (J/K/mol), (K)
    constant_expressions  = '0.64  8.3144598  1073.0'
    function = 'D_O2/R/T*101325'
  [../]

  [./vacancyDiffFluxCoefficient]
    type  = ParsedMaterial
    block = 'PT_LSM_TET4 PT_TPB_TET4'
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    f_name = 'vacDiffFluxCoef'
    constant_names        = 'D_O      a         NA'     # (cm^2/s), (cm), (1/mol), .
    constant_expressions  = '7.5e-7   3.893e-8  6.022e23'
    function = 'D_O/(a^3)/NA*1e6'
  [../]

  [./vacancyDriftFluxCoefficient]
    type  = GenericConstantMaterial
    block = 'PT_YSZ_TET4 PT_TPB_TET4 PT_IYSZ_TET4'
    prop_names  = 'sigma_YSZ' # (S/cm)
    prop_values = '4e-2'
  [../]

  [./holeConductivity]
    type  = GenericConstantMaterial
    block = 'PT_LSM_TET4 PT_TPB_TET4' # when you use previously generated mesh, you need to use LSM block.
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    prop_names  = 'sigma_e' # (S/cm)
    prop_values = '1.43e5' # Conductivity at 293 K
  [../]

  [./CounterElectrode]
    type  = ParsedMaterial
    block = 'PT_YSZ_TET4 PT_TPB_TET4 PT_IYSZ_TET4'
    # block = 'PT_YSZ_C_TET4 PT_YSZ_A_TET4 PT_TPB_C_TET4 PT_TPB_A_TET4' # Full cell
    f_name = 'eff_pO2_CE'
    # constant_names        = 'R          T        z     F'        # (J/K/mol), (K),  (# of electron), (C/mol)
    # constant_expressions  = '8.3144598  1073.0   4   96485.3329'
    # postprocessor_names = 'numerical_phi_YSZ'
    # in terms of effective Oxygen partial pressure (atm)
    # function = 'if(numerical_phi_YSZ=0, exp(z * F / R / T * (-1.07277)), exp(z * F / R / T * (-(1.03665 + numerical_phi_YSZ))))'
    function = '7e-21' # (atm)
    # Make values see at output(*.e) file
    # outputs = exodus
  [../]

  [./tpbActivity]
    type  = ParsedMaterial
    block = 'PT_TPB_TET4'
    f_name = 'tpbActivity_S0'
    # constant_names        = 'R          T        z     F          x_CrO2OH2   x_H2O  V_molar_Cr2O3'
    #                     # (J/K/mol),  (K), (# of e-), (C/mol),    (non),      (non), (cm^3)
    # constant_expressions  = '8.3144598  1073.0   6   96485.3329   2.5e-9       0.01   29.12'
    function = 's0_tpb_init'
                # Here, "CrDeposit * time" can be changed depending on how much do I want to make a Cr deposition.
                # This if-else statement is needed to check if tpbACtivity is less than zero.
    args = 's0_tpb_init'
    # material_property_names = 'VolCr time'
    # postprocessor_names = ''
    # Make values see at output(*.e) file
    # outputs = exodus
  [../]

  # Provide various time stepping quantities as material properties (dt / t / time_step)
  # Need Phase_field Module, so if this does not work, then "Check your Makefile"
  [./ProvidingTime]
    type = TimeStepMaterial
    block = 'PT_TPB_TET4'
  [../]
[]

#==========================================================================#

[Variables]
  [./p_O2]
    block = 'PT_PORE_TET4 PT_TPB_TET4'
    initial_condition = 0.21 # (atm)
    scaling = 1e4
  [../]

  [./V_O]
    block = 'PT_LSM_TET4 PT_TPB_TET4'
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    initial_condition = 2.580947226225166e-08 # (.) pO2 = 0.21 atm
    scaling = 1e8
  [../]

  [./phi_LSM]
    block = 'PT_LSM_TET4 PT_TPB_TET4' # when you use previously generated mesh, you need to use LSM block.
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    initial_condition = 1.07272 # (V) Just for condition: Each subdomain must contain at least one Kernel.
    # initial_condition = 1.02845 # (V) Just for condition: Each subdomain must contain at least one Kernel.
    scaling = 1e2 # 1e2 for cylinder
  [../]

  [./phi_YSZ]
    block = 'PT_YSZ_TET4 PT_TPB_TET4 PT_IYSZ_TET4'
    # initial_condition = 0.00000 # (V)
    initial_condition = 0.03607 # (V)
    scaling = 1e4
  [../]
[]

#==========================================================================#

# [UserObjects]
#   [./CheckPSmode]
#     type = Terminator
#     expression = '(abs(j_YSZ_bottom) < 1e-2)' # OR condition
#     fail_mode = HARD
#     error_level = INFO
#     message = 'Current density dropped below limit : 1e-3 or Remainder of an integer division for time step is not 0'
#     execute_on = timestep_end
#   [../]
# []

#==========================================================================#

[Kernels]
  [./gasDiffusion]
    type  = DiffMatKernel
    block = 'PT_PORE_TET4 PT_TPB_TET4'
    variable  = p_O2
    diff_coef = 'gasDiffFluxCoef'
  [../]

  [./vacancyDiffusion]
    type  = DiffMatKernel
    block = 'PT_LSM_TET4 PT_TPB_TET4'
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    variable  = V_O
    diff_coef = 'vacDiffFluxCoef'
  [../]

  [./holeDrift]
    type  = DiffMatKernel
    block = 'PT_LSM_TET4 PT_TPB_TET4' # when you use previously generated mesh, you need to use LSM block.
    # block = 'PT_NI_TET4 PT_TPB_TET4'
    variable  = phi_LSM
    diff_coef = 'sigma_e'
  [../]

  [./vacancyIonicDrift]
    type  = DiffMatKernel
    block = 'PT_YSZ_TET4 PT_TPB_TET4 PT_IYSZ_TET4'
    variable  = phi_YSZ
    diff_coef = 'sigma_YSZ'
  [../]

  [./tpbReactionOxygenPore]
    type = CoupledTPBOxygenPressurePoreQS
    block = 'PT_TPB_TET4'
    variable = p_O2
    phi_YSZ = phi_YSZ
    phi_LSM = phi_LSM
    pO2_CE = 'eff_pO2_CE'
    s0     = 'tpbActivity_S0'
  [../]

  [./tpbReactionPotentialYSZ]
    type = CoupledTPBPotentialYSZQS
    block = 'PT_TPB_TET4'
    variable = phi_YSZ
    p_O2 = p_O2
    phi_LSM = phi_LSM
    pO2_CE = 'eff_pO2_CE'
    s0     = 'tpbActivity_S0'
  [../]

  [./tpbReactionPotentialLSM]
    type = CoupledTPBPotentialLSMQS
    block = 'PT_TPB_TET4'
    variable = phi_LSM
    p_O2 = p_O2
    phi_YSZ = phi_YSZ
    pO2_CE = 'eff_pO2_CE'
    s0     = 'tpbActivity_S0'
  [../]
[]

#==========================================================================#

[InterfaceKernels]
  [./interfaceSurfaceExchangeFullyCoupled]
    type = InterfaceSurfExchangeFullyCoupled
    variable = p_O2
    neighbor_var = V_O
    boundary = 'interface_12'
    k = 6.14e-6 # (cm/s)
  [../]

  [./interfaceChargeTransferFullyCoupled]
    type = InterfaceChargeTransferFullyCoupledQS
    variable = V_O
    neighbor_var = phi_YSZ
    boundary = 'interface_23 interface_26'
    j0 = 0.193  # (A/cm^2)
    function_phi_LSM = 'funcPotentialLSM'
  [../]
[]

#==========================================================================#

[BCs]
  [./oxygenPartialPressure_top]
    type = DirichletBC
    variable = p_O2
    boundary = 'SF_PORE_WITH_XMIN'
    # boundary = 'SF_PORE_WITH_ZMAX' # cylinder geometry
    value = 0.21 # (atm)
  [../]

  [./potentialLSM_top]
    type = FunctionDirichletBC
    variable = phi_LSM
    boundary = 'SF_LSM_WITH_XMIN'
    # boundary = 'SF_NI_WITH_ZMAX' # cylinder geometry
    function = 'funcPotentialLSM' # (V)
  [../]

  [./potentialYSZ_bottom]
    type = DirichletBC
    variable = phi_YSZ
    boundary = 'SF_YSZ_WITH_XMAX'
    # boundary = 'SF_YSZ_WITH_ZMIN' # cylinder geometry
    # value = 0.00000 # (V)
    value = 0.03607 # (V)
  [../]
[]

#==========================================================================#

[AuxVariables]
  [./aux_eta_tpb]
    block = 'PT_TPB_TET4'
  [../]

  [./s0_tpb_init]
    block = 'PT_TPB_TET4'
    initial_condition = 5318
  [../]
[]

[AuxKernels]
  [./eta_tpb]
    type = ParsedAux
    variable = aux_eta_tpb
    block = 'PT_TPB_TET4'
    function = '-R * T / 4 / F * log(7e-21 / p_O2) - (phi_LSM - phi_YSZ)'
    constant_names = 'R T F'
    constant_expressions = '8.3144598 1073 96485.33289' # (J/K/mol), (K), (C/mol)
    args = 'p_O2 phi_LSM phi_YSZ'
  [../]
[]

#==========================================================================#

[Postprocessors]
  [./I_YSZ_bottom]
    # type = SideFluxIntegral
    type = SideDiffusiveFluxIntegral
    variable = phi_YSZ
    diffusivity = 'sigma_YSZ'
    boundary = 'SF_YSZ_WITH_XMAX'
    # boundary = 'SF_YSZ_WITH_ZMIN' # cylinder geometry
    outputs = 'console csv'
  [../]

  [./j_YSZ_bottom]
    # type = SideFluxAverage
    type = SideDiffusiveFluxAverage
    variable = phi_YSZ
    diffusivity = 'sigma_YSZ'
    boundary = 'SF_YSZ_WITH_XMAX'
    # boundary = 'SF_YSZ_WITH_ZMIN' # cylinder geometry
    outputs = 'console csv'
  [../]

  [./numerical_phi_YSZ]
    type = AverageNodalVariableValue
    # block = 'PT_YSZ_TET4' # Combined YSZ in Full cell
    boundary = 'SF_YSZ_WITH_XMAX'
    # boundary = 'SF_YSZ_WITH_ZMIN' # cylinder geometry, Half cell simulation
    # boundary = 'SF_CounterElectrode' # Separated YSZ in Full cell simulation
    variable = phi_YSZ
    execute_on = 'initial timestep_begin nonlinear timestep_end' # Need this line, because Postprocessor is not being executed at the 0th timestep.
    # execute_on = 'initial timestep_end' # Need this line, because Postprocessor is not being executed at the 0th timestep.
  [../]

  [./phi_LSM_model]
    type = FunctionValuePostprocessor
    function = 'funcPotentialLSM'
    outputs = 'console csv'
    execute_on = 'initial timestep_begin nonlinear timestep_end' # Need this line, because Postprocessor is not being executed at the 0th timestep.
  [../]

  [./eta_total]
    type = FunctionValuePostprocessor
    function = 'funcOverpotential'
    outputs = 'console csv'
    execute_on = 'initial timestep_begin nonlinear timestep_end' # Need this line, because Postprocessor is not being executed at the 0th timestep.
  [../]

  [./dt] # This is needed for Terminator condition
    type = TimestepSize
  [../]
[]

#==========================================================================#

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    #off_diag_column = p_O2
    #off_diag_row = phi_YSZ
    petsc_options = '-snes_converged_reason -ksp_converged_reason'
    #petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type -snes_type'
    #petsc_options_value = 'lu       1e-6                 ds           test'
    #petsc_options_iname = '-snes_type'
    #petsc_options_value = 'test'
  [../]
[]

#==========================================================================#

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 1.0
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-8
  nl_abs_tol = 1e-6
  # l_tol = 1e-04
  # l_abs_step_tol = -1
  l_max_its = 1000

  [./TimeStepper]
    type = FunctionDT
    function = funcTimeStepper
  [../]

  # solve_type = 'NEWTON'
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = 'hypre boomeramg 50'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type'
  # petsc_options_value = '101 bjacobi'

  solve_type = 'NEWTON'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  petsc_options_value = '101 bjacobi ilu'
[]

#==========================================================================#

[Outputs]
  exodus = true
  file_base = outputs/cathode_intial
  append_date = true
  append_date_format = '%Y-%m-%d-%R-%S'
  # execute_on = 'initial final'
  [csv]
    type = CSV
    file_base = outputs/cathode_test
    append_date = true
    append_date_format = '%Y-%m-%d-%R-%S'
  []
  # [exodus1]
  #   type = Exodus
  #   file_base = outputs/cathode_test1
  #   append_date = true
  #   append_date_format = '%Y-%m-%d-%R-%S'
  #   sync_only = true
  #   sync_times = '0 200 400 600 800 1000'
  #   output_material_properties = true
  #   show_material_properties = 'tpbActivity_S0 VolCr'
  # []
  # [exodus2]
  #   type = Exodus
  #   file_base = outputs/cathode_test2
  #   append_date = true
  #   append_date_format = '%Y-%m-%d-%R-%S'
  #   output_material_properties = true
  #   show_material_properties = 'tpbActivity_S0 VolCr'
  #   execute_on = 'initial final'
  # []
[]

[Debug]
  show_var_residual_norms = true
[]
