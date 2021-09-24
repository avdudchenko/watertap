import sys
import os
import time

from proteuslib.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep


def append_costing_outputs(m, outputs, units_to_cost):
    for unit_name in units_to_cost:
        for cost_type in ['capital_cost', 'operating_cost']:
            unit = getattr(m.fs, unit_name)
            cost = getattr(unit.costing, cost_type)

            outputs['%s_%s' % (cost_type, unit_name)] = cost

    return outputs


def run_analysis(case_num, nx, RO_type):

    desal_kwargs = {'has_desal_feed': False, 'is_twostage': True, 'has_ERD': True,
                    'RO_type': RO_type, 'RO_base': 'TDS', 'RO_level': 'detailed'}

    if case_num == 1:
        # ================================================================
        # Single Stage RO
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_single_stage as fs_single_stage

        desal_kwargs['is_twostage'] = False

        m = fs_single_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        # sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.95, nx)
        sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW
        outputs['Saturation Index'] = m.fs.desal_saturation.saturation_index
        outputs['Pump Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
        outputs = append_costing_outputs(m, outputs, ['RO', 'pump_RO', 'ERD'])

        output_filename = 'output/fs_single_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_single_stage.optimize

    elif case_num == 2:
        # ================================================================
        # Two Stage RO
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as fs_two_stage 

        m = fs_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        # sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.95, nx)
        sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.3, 0.95, nx)
        sweep_params['Max Pressure'] = LinearSample(m.fs.max_allowable_pressure, 300e5, 75e5, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW
        outputs['Saturation Index'] = m.fs.desal_saturation.saturation_index
        outputs['RO-1 Pump Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
        outputs['RO-2 Pump Pressure'] = m.fs.pump_RO2.control_volume.properties_out[0].pressure
        outputs = append_costing_outputs(m, outputs, ['RO', 'pump_RO', 'RO2', 'pump_RO2', 'ERD'])

        output_filename = 'output/fs_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_two_stage.optimize

    elif case_num == 3:
        # ================================================================
        # NF No Bypass
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_no_bypass as fs_NF_no_bypass

        m = fs_NF_no_bypass.solve_flowsheet()

        sweep_params = {}
        sweep_params['NF Recovery'] = LinearSample(m.fs.NF.recovery_mass_phase_comp[0, 'Liq', 'H2O'], .3, .95, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW
        outputs['Saturation Index'] = m.fs.pretrt_saturation.saturation_index

        output_filename = 'output/fs_NF_no_bypass/results_%d.csv' % case_num

        opt_function = fs_NF_no_bypass.simulate

    elif case_num == 4:
        # ================================================================
        # NF Two Stage
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_two_stage as fs_NF_two_stage

        m = fs_NF_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], .01, .99, nx)

        outputs = {}
        outputs['Ca Removal'] = m.fs.removal_Ca

        output_filename = 'output/fs_NF_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_NF_two_stage.optimize

    elif case_num == 5:
        # ================================================================
        # NF Two Stage
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_two_stage as fs_NF_two_stage

        m = fs_NF_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], 0.25, 1.0, 4)
        # sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)
        sweep_params['RO Recovery'] = LinearSample(m.fs.RO_recovery, 0.5, 0.95, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW

        output_filename = 'output/fs_NF_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_NF_two_stage.optimize

    elif case_num == 6:
        # ================================================================
        # NF Two Stage
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_two_stage as fs_NF_two_stage

        m = fs_NF_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

        outputs = {}
        outputs['NF Split Fraction'] = m.fs.splitter.split_fraction[0, 'bypass']

        output_filename = 'output/fs_NF_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_NF_two_stage.optimize

    elif case_num == 7:
        # ================================================================
        # NF Two Stage
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF_two_stage as fs_NF_two_stage

        m = fs_NF_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        sweep_params['NF Split Fraction'] = LinearSample(m.fs.splitter.split_fraction[0, 'bypass'], 0.25, 1.0, 4)
        sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW

        output_filename = 'output/fs_NF_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_NF_two_stage.optimize

    elif case_num == 8:
        # ================================================================
        # Softening
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening as fs_softening

        m = fs_softening.solve_flowsheet()

        sweep_params = {}
        sweep_params['Lime Dosing'] = LinearSample(m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol, 0.0, 0.128, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW
        outputs['Ca Removal'] = m.fs.removal_Ca
        outputs['Mg Removal'] = m.fs.removal_Mg

        output_filename = 'output/fs_softening/results_%d.csv' % case_num

        opt_function = fs_softening.simulate

    elif case_num == 9:
        # ================================================================
        # Softening Two Stage
        # ================================================================
        import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening_two_stage as fs_softening_two_stage

        m = fs_softening_two_stage.optimize_flowsheet(**desal_kwargs)

        sweep_params = {}
        sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, nx)

        outputs = {}
        outputs['LCOW'] = m.fs.costing.LCOW
        outputs['Lime Dosing'] = m.fs.stoich_softening_mixer_unit.lime_stream_state[0].flow_mol

        output_filename = 'output/fs_softening_two_stage/results_%d_%sRO.csv' % (case_num, desal_kwargs['RO_type'])

        opt_function = fs_softening_two_stage.optimize

    else:
        raise ValueError('case_num = %d not recognized.' % (case_num))


    global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
                                     optimize_function=opt_function,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local')#,
                                     #interpolate_nan_outputs=True) FIXME: This only works with 2D sweeps

    return global_results, sweep_params


if __name__ == "__main__":

    # Start MPI communicator
    comm, rank, num_procs = _init_mpi()

    # Get the case number to run
    try:
        case_num = int(sys.argv[1])
    except:
        # Default to running case 1
        case_num = 1

    # Get the default number of discretization points
    try:
        nx = int(sys.argv[2])
    except:
        # Default to a 4-point discretization
        nx = 4

    # Get the default RO unit to use (1D or 0D)
    try:
        RO_type = sys.argv[3]
    except:
        # Default to RO_type='0D'
        RO_type = '0D'

    tic = time.time()
    global_results, sweep_params = run_analysis(case_num, nx, RO_type)
    print(global_results)
    toc = time.time()

    if rank == 0:
        total_samples = 1

        for k, v in sweep_params.items():
            total_samples *= v.num_samples
            
        print('Finished case_num = %d.' % (case_num))
        print('Processed %d swept parameters comprising %d total points.' % (len(sweep_params), total_samples))
        print('Elapsed time = %.1f s.' % (toc-tic))
