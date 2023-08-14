###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from idaes.core.solvers import get_solver

from watertap.tools.parameter_sweep.parameter_sweep import (
    ParameterSweep,
    RecursiveParameterSweep,
)

from watertap.tools.parameter_sweep.parameter_sweep_differential import (
    DifferentialParameterSweep,
)

from watertap.tools.parameter_sweep.parameter_sweep_reader import ParameterSweepReader

from watertap.tools.parameter_sweep.parameter_sweep_differential import (
    DifferentialParameterSweep,
)

from watertap.tools.analysis_tools.loop_tool.data_merging_tool import *

from watertap.tools.parallel.parallel_manager_factory import (
    has_mpi_peer_processes,
    get_mpi_comm_process,
)
import copy
import os
import h5py
import numpy as np

import watertap.tools.MPI as MPI


import multiprocessing


def get_working_dir():
    cwd = os.getcwd()
    return cwd


class loopTool:
    def __init__(
        self,
        loop_file,
        solver=None,
        build_function=None,
        initialize_function=None,
        optimize_function=None,
        number_of_subprocesses=None,
        probe_function=None,
        save_name=None,
        saving_dir=None,
        execute_simulations=True,
        custom_do_param_sweep=None,
        custom_do_param_sweep_kwargs=None,
        h5_backup=None,
    ):
        """
        Loop tool class that runs iterative paramter sweeps

        Arguments:
            loop_file : .yaml config file that contains iterative loops to run
            solver : solver to use in model, default uses watertap solver
            build_function : function to build unit model
            initialize_function : function for intilization of th eunit model
            optimize_function : function for solving model
            probe_function : Function to probe if a solution should be attempted or not
            save_name : name to use when saving the file with
            save_dir : directory to save the file in
            number_of_subprocesses : user defined number of subprocesses to use for parallel run, defaults to either
            max number of logical cores, if set to False, will disable MPI and set number_of_subpressess to 0
            custom_do_param_sweep : custom param function (refer to parameter sweep tool)
            custom_do_param_sweep_kwargs : custom parm kwargs (refer to parameter sweep tool)
            execute_simulations : of looptool should execute simulations upon setup,
                                    other user can call build_run_dict, and run_simulations call manually
            h5_backup : Set location for back up file, if set to False, no backup will be created, otherwise backup will be autocreated

        """

        self.loop_file = loop_file
        self.solver = solver

        self.build_function = build_function
        self.initialize_function = initialize_function
        self.optimize_function = optimize_function
        self.probe_function = probe_function

        self.save_name = save_name
        self.data_dir = saving_dir

        self.number_of_subprocesses = number_of_subprocesses

        self.setup_multi_processing()

        self.test_mode = False

        self.custom_do_param_sweep = custom_do_param_sweep
        self.custom_do_param_sweep_kwargs = custom_do_param_sweep_kwargs
        self.h5_backup_location = h5_backup
        if execute_simulations:
            self.build_run_dict()
            self.run_simulations()

    def build_run_dict(self, test_setups=False):
        """
        This builds the dict that will be used for simulatiuons

        Arguments:
            test_setups : test if configuraiton will intilaize, but not run the simulatuons
            create_directories : if save directories should be created on the fly
        """

        loop_dict = ParameterSweepReader()._yaml_to_dict(self.loop_file)
        self.sweep_directory = {}
        for key, loop in loop_dict.items():
            self.sweep_directory[key] = {}
            self.init_sim_options()
            self.save_dir = self.save_dir + "/" + key
            self.h5_directory = key
            loop_type = self.get_loop_type(loop)
            self.options = loop
            self.sweep_directory[key], dir = self.build_sweep_directories(
                loop[loop_type],
                loop_type,
                self.sweep_directory[key],
                self.save_dir,
                self.h5_directory,
            )
            if test_setups:
                self.execute_sweep(self.sweep_directory[key], False)

    def run_simulations(self):
        """runs the simulations crated in build_run_dict"""

        self.execute_sweep(self.sweep_directory)

    def run_model_test(self):
        """function to run a single test, with out param sweep"""
        self.test_mode = True
        self.execute_sweep(self.sweep_directory)
        self.test_mode = False

    def get_loop_type(self, loop):
        """containst types of loop options that are extracted from yaml file
        build_loop - options that are passed into build function
        init_loop - options that are passed into initialization function
        optimize_loop - options that are passed into optimize function
        sweep_param_loop - this will creates the paramters that should be sweeped over
        diff_param_loop - parmaters for differential sweeps
        """

        if "build_loop" in loop:
            loop_type = "build_loop"
        elif "init_loop" in loop:
            loop_type = "init_loop"
        elif "optimize_loop" in loop:
            loop_type = "optimize_loop"
        elif "sweep_param_loop" in loop:
            loop_type = "sweep_param_loop"
        elif "diff_param_loop" in loop:
            loop_type = "diff_param_loop"
        else:
            loop_type = None
        return loop_type

    def get_loop_key(self, loop, loop_type):
        """function for finding loop type"""
        for key in loop.keys():
            if key != loop_type:
                return key

    def build_sweep_directories(
        self, loop, loop_type, sweep_directory, cur_dir, cur_h5_dir
    ):
        """this creats the loop directory dict, which is then used to run
        the paramter sweep"""
        if loop_type != None:
            loop_type_recursive = self.get_loop_type(loop)
            loop_key_current = self.get_loop_key(loop, loop_type)
            if loop_type_recursive != None:
                sweep_directory[loop_key_current] = {}
                for loop_value in loop[loop_key_current]:
                    if isinstance(loop_value, dict):
                        loop_value = str(loop_value).strip("{}")
                    sweep_directory[loop_key_current][loop_value] = {}
                    cur_dir = self.update_dir_path(
                        cur_dir, loop_key_current, loop_value
                    )
                    cur_h5_dir = self.update_dir_path(
                        cur_h5_dir, loop_key_current, loop_value
                    )
                    self.update_sim_options(
                        loop_type, loop_key_current, loop_value, loop
                    )
                    (
                        sweep_directory[loop_key_current][loop_value],
                        cur_dir,
                    ) = self.build_sweep_directories(
                        loop[loop_type_recursive],
                        loop_type_recursive,
                        sweep_directory[loop_key_current][loop_value],
                        cur_dir,
                        cur_h5_dir,
                    )
            else:
                for loop_value in loop:
                    local_dir = self.update_dir_path(cur_dir, loop_value, value=None)
                    h5_local_dir = self.update_dir_path(
                        cur_h5_dir, loop_value, value=None
                    )
                    self.update_sim_options(loop_type, loop_value, loop, None)
                    # creates directory dict with all options
                    if bool(self.sweep_params):
                        sweep_directory[loop_value] = {
                            "simulation_setup": {
                                "dir": local_dir,
                                "h5dir": h5_local_dir,
                                "build_defaults": copy.deepcopy(self.build_defaults),
                                "init_defaults": copy.deepcopy(self.init_defaults),
                                "optimize_defaults": copy.deepcopy(
                                    self.optimize_defaults
                                ),
                                "sweep_params": copy.deepcopy(self.sweep_params),
                                "num_samples": copy.deepcopy(self.num_samples),
                                "expected_num_samples": copy.deepcopy(
                                    self.expected_num_samples
                                ),
                                "force_rerun": copy.deepcopy(self.force_rerun),
                                "min_num_samples": copy.deepcopy(self.min_num_samples),
                                "diff_params": copy.deepcopy(self.diff_params),
                                "diff_samples": copy.deepcopy(self.diff_samples),
                                "original_options_dict": copy.deepcopy(self.options),
                            }
                        }
                    else:
                        print("loop empty did not add", loop_value)
        return sweep_directory, cur_dir

    def update_dir_path(self, cur_dir, key, value):
        """creates directory pathing for h5 file strucutre"""
        if value == None:
            cur_dir = cur_dir + "/" + str(key)
        else:
            if cur_dir.find(key) > 0:
                cur_dir = cur_dir[: (cur_dir.find(key) - 1)]
            cur_dir = cur_dir + "/" + str(key)
            cur_dir = cur_dir + "/" + str(value)
        return cur_dir

    def update_sim_options(self, loop_type, loop_key, loop_value, loop):
        """used to update simulation options for each specific loop"""

        if loop_type == "build_loop":
            self.build_defaults.update(self.get_loop_params(loop_key, loop_value, loop))
        elif loop_type == "init_loop":
            self.init_defaults.update(self.get_loop_params(loop_key, loop_value, loop))
        elif loop_type == "optimize_loop":
            self.optimize_defaults.update(
                self.get_loop_params(loop_key, loop_value, loop)
            )
        elif loop_type == "sweep_param_loop":
            self.sweep_params = {}
            (
                self.sweep_params,
                self.num_samples,
                self.expected_num_samples,
                self.min_num_samples,
                self.force_rerun,
            ) = self.get_sweep_params(loop_key, loop_value)

        elif loop_type == "diff_param_loop":
            self.sweep_params = {}
            self.diff_params = {}
            if loop_key != "sweep_reference_params":
                (
                    self.sweep_params,
                    self.num_samples,
                    self.diff_params,
                    self.diff_samples,
                    self.expected_num_samples,
                    self.min_num_samples,
                    self.force_rerun,
                ) = self.get_diff_params(loop_key, loop_value)

    def get_loop_params(self, key, loop_value, loop):
        if "case" in key:
            return loop[key][loop_value]
        elif ":" in str(loop_value):
            lp = loop_value.split(":")
            return {key: {eval(lp[0]): float(lp[1])}}
        else:
            return {key: loop_value}

    def setup_multi_processing(self):
        if self.number_of_subprocesses == False:
            self.mpi_comm = False
            self.number_of_subprocesses = 1
        else:
            if has_mpi_peer_processes():
                self.mpi_comm = get_mpi_comm_process()
            else:
                self.mpi_comm = False

    def get_diff_params(self, key, loop_value):
        """creates dict for differntial sweep
        -creates dict for diff_spec, these are differnetial samples that would be sampled for each
        reference design
        -creaste sweep_reference-dict which will contain paramters that will be sweeped over to generate reference
        values
        """
        sweep_params = {}
        diff_params = {}
        sweep_samples = 1
        min_num_samples = None
        force_rerun = None
        # print(loop_value, key)
        if "diff_mode" in loop_value[key]:
            diff_samples = loop_value[key]["num_samples"]
            diff_params[loop_value[key]["param"]] = loop_value[key]
            expected_num_samples = loop_value[key]["num_samples"]
            min_num_samples = loop_value[key].get("min_num_samples")
            force_rerun = loop_value[key].get("rerun")
        else:
            for key, values in loop_value[key].items():
                # print(values)
                diff_samples = values["num_samples"]
                diff_params[values["param"]] = values
                expected_num_samples = values["num_samples"]
                min_num_samples = loop_value[key].get("min_num_samples")
                force_rerun = loop_value[key].get("rerun")
        sweep_samples = loop_value["sweep_reference_params"]["num_samples"]
        for key, values in loop_value["sweep_reference_params"].items():
            if key != "num_samples" and key != "min_num_samples":
                sweep_params[values["param"]] = values
                if "num_samples" not in values:
                    sweep_params[values["param"]]["num_samples"] = sweep_samples
        return (
            sweep_params,
            sweep_samples,
            diff_params,
            diff_samples,
            expected_num_samples,
            min_num_samples,
            force_rerun,
        )

    def get_sweep_params(self, key, loop_value):
        if "type" in loop_value[key]:
            num_samples = loop_value[key]["num_samples"]
            try:
                expected_num_samples = loop_value[key]["expected_num_samples"]
            except KeyError:
                expected_num_samples = num_samples
            # try:
            min_num_samples = loop_value[key].get("min_num_samples")
            force_rerun = loop_value[key].get("rerun")
            param = loop_value[key]["param"]
            # except:
            return (
                {param: loop_value[key]},
                num_samples,
                expected_num_samples,
                min_num_samples,
                force_rerun,
            )
        else:
            params = {}
            num_samples = 1
            expected_num_samples = None
            min_num_samples = None
            force_rerun = None
            for key, values in loop_value[key].items():
                if key == "expected_num_samples":
                    expected_num_samples = values
                elif key == "min_num_samples":
                    min_num_samples = values
                elif key == "rerun":
                    force_rerun = values
                else:
                    param = values["param"]
                    params[param] = values
                    num_samples = num_samples * values["num_samples"]
            if expected_num_samples == None:
                expected_num_samples = num_samples
            return (
                params,
                num_samples,
                expected_num_samples,
                min_num_samples,
                force_rerun,
            )

    def init_sim_options(self):
        """resets simulations options"""

        self.build_defaults = {}
        self.init_defaults = {}
        self.optimize_defaults = {}
        self.outputs = None
        self.sweep_params = {}
        self.diff_params = {}
        self.diff_samples = 0

        self._create_save_directory(self.data_dir)
        self.save_dir = self.data_dir + "/output"

        self._create_save_directory(self.save_dir)
        self.save_dir = self.save_dir + "/" + self.save_name
        self.h5_file_location_default = self.save_dir
        self.h5_directory = ""

    def execute_sweep(
        self,
        sweep_directory,
    ):
        """runs through sweep directory and execute each simulation
        unless in test mode, in which only try build, init, and solve"""
        for key, value in sweep_directory.items():
            if key != "simulation_setup":
                self.execute_sweep(value)
            elif self.test_mode:
                self.test_setup(value)
                break
            else:
                self.execute_param_sweep_run(value)

    def execute_param_sweep_run(self, value):
        """this executes the parameter sweep
        if no data exists in the back upfile, or user forces exceution
        this defined by check solution exists function
        """
        self.setup_param_sweep(value)
        solution_test = self.check_solution_exists()
        if solution_test:
            self.build_sim_kwargs()
            self.sweep_params = value["sweep_params"]
            if value["diff_params"] == {}:
                self.run_class_parameter_sweep()
            else:
                self.differential_sweep_specs = value["diff_params"]
                self.diff_samples = value["diff_samples"]
                self.run_diff_parameter_sweep()

    def setup_param_sweep(self, value):
        """set up variables before a sweep run with parmater sweep
        tool, resets any of prior options"""
        self.init_sim_options()
        self.options = value["original_options_dict"]
        self.build_default = value["build_defaults"]
        self.outputs = self.options.get("outputs")
        self.optimize_defaults = value["optimize_defaults"]

        self.init_defaults = value["init_defaults"]
        self.save_dir = value["dir"]
        self.h5_directory = value["h5dir"]
        val = self.h5_directory.split("/")[0]

        self.h5_file_location = (
            self.h5_file_location_default + "_analysisType_" + str(val) + ".h5"
        )
        # resets it if file name changes
        if (
            self.h5_backup_location is not None
            and self.h5_file_location not in self.h5_backup_location
        ):
            self.h5_backup_location = None
        self.num_samples = value["num_samples"]
        self.expected_num_samples = value["expected_num_samples"]
        self.force_rerun = value["force_rerun"]
        self.min_num_samples = value["min_num_samples"]

    def build_sim_kwargs(self):
        """here we build all the kwargs option
        first loading defaults provided by default calls, followed
        by those in the loops, overriding any defaults or adding new
        options to kwargs"""
        # check if user wants to reinit befoure sweep.
        if "reinitialize_before_sweep" in self.options:
            self.reinitialize_before_sweep = self.options["reinitialize_before_sweep"]
        else:
            self.reinitialize_before_sweep = False
        # generated combined build kwargs (default + loop)
        self.combined_build_defaults = {}  # self.build_default
        self.combined_build_defaults.update(self.options.get("build_defaults", {}))
        self.combined_build_defaults.update(self.build_default)
        # generated combined optimize kwargs (default + loop)
        self.combined_optimize_defaults = {}
        self.combined_optimize_defaults.update(
            self.options.get("optimize_defaults", {})
        )
        self.combined_optimize_defaults.update(self.optimize_defaults)
        # generated combined init  kwargs (default + loop)
        self.combined_init_defaults = {}
        self.combined_init_defaults.update(self.options.get("init_defaults", {}))
        self.combined_init_defaults.update(self.init_defaults)

    def test_setup(self, value):
        """test if passed in kawrgs will work
        ToDO: add parmas in paramters sweep to test"""

        self.setup_param_sweep(value)
        self.build_sim_kwargs()
        solver = get_solver()
        reinit_kwarg = {"solver": solver}
        reinit_kwarg.update(self.combined_init_defaults)
        m = self.build(**self.combined_build_defaults)
        self.initialize_system(m, **reinit_kwarg)
        self.solve(m, solver=solver)

    def _check_solution_exists(self):
        """hidden function to check if solution
        exists, and create and h5 file for data storage"""
        self.h5_backup_location = create_backup_file(
            self.h5_file_location, self.h5_backup_location, self.h5_directory
        )
        create_h5_file(self.h5_file_location)

        sucess_feasible = self.expected_num_samples
        if self.force_rerun == True:
            print("Forceing rerun", self.force_rerun)
            run_sweep = True
        else:
            run_sweep = merge_data_into_file(
                self.h5_file_location,
                self.h5_backup_location,
                self.h5_directory,
                expected_solved_values=self.expected_num_samples,
                min_solve_values=self.min_num_samples,
                force_rerun=self.force_rerun,
            )
        return run_sweep

    def check_solution_exists(self):
        """check if solution exists, ensuring to use only
        rank 0 if user running MPI, otherwise do direct
        solution check"""
        print("Checking if solution exists", self.h5_file_location, self.h5_directory)
        self.cur_h5_file = (self.h5_file_location, self.h5_directory)
        if self.mpi_comm != False:
            self.mpi_comm.Barrier()
            results = np.empty(self.mpi_comm.Get_size(), dtype=bool)

            results[:] = True
            if self.mpi_comm.Get_rank() == 0:
                success = self._check_solution_exists()
                results[:] = success
                print("Got test_result", results)
            self.mpi_comm.Bcast(results, root=0)
            self.mpi_comm.Barrier()

            success = results[self.mpi_comm.Get_rank()]
        else:
            success = self._check_solution_exists()
        return success

    def _create_save_directory(self, save_dir):
        """used to create a save directory (outputs)"""
        try:
            os.mkdir(save_dir)
        except OSError as error:
            pass

    def run_class_parameter_sweep(self):
        """setup and run paramer sweep"""

        # get solver if not provided by user
        if self.solver is None:
            solver = get_solver()
        else:
            solver = self.solver

        # add solver to init kwargs
        self.combined_init_defaults.update({"solver": solver})
        # add solver to optimize kwarg
        self.combined_optimize_defaults.update({"solver": solver})

        # setup parmater sweep tool
        ps_kwargs = {}
        ps_kwargs["csv_results_file_name"] = None
        ps_kwargs["h5_results_file_name"] = self.cur_h5_file[0]
        ps_kwargs["h5_parent_group_name"] = self.cur_h5_file[1]

        ps_kwargs["optimize_function"] = self.optimize_function
        ps_kwargs["optimize_kwargs"] = self.combined_optimize_defaults

        ps_kwargs["reinitialize_function"] = self.initialize_function
        ps_kwargs["reinitialize_kwargs"] = self.combined_init_defaults
        ps_kwargs["reinitialize_before_sweep"] = self.reinitialize_before_sweep

        ps_kwargs["custom_do_param_sweep"] = self.custom_do_param_sweep
        ps_kwargs["custom_do_param_sweep_kwargs"] = self.custom_do_param_sweep_kwargs

        ps_kwargs["probe_function"] = self.probe_function

        ps_kwargs["number_of_subprocesses"] = self.number_of_subprocesses

        ps = ParameterSweep(**ps_kwargs)
        ps.parameter_sweep(
            self.build_function,
            ParameterSweepReader()._dict_to_params,
            build_outputs=None,
            num_samples=self.num_samples,
            build_model_kwargs=self.combined_build_defaults,
            build_sweep_params_kwargs={"input_dict": self.sweep_params},
        )

    def run_diff_parameter_sweep(self):
        """setup and run diff paramer sweep
        fix once the diff paramter tool is updated to new version."""

        if self.solver is None:
            solver = get_solver()
        else:
            solver = self.solver

        # add solver to init kwargs
        self.combined_init_defaults.update({"solver": solver})
        # add solver to optimize kwarg
        self.combined_optimize_defaults.update({"solver": solver})

        # legacy, need to build m
        m = self.build_function(**self.combined_build_defaults)

        # setup parmater sweep tool
        ps_kwargs = {}
        ps_kwargs["csv_results_file_name"] = None
        ps_kwargs["h5_results_file_name"] = self.cur_h5_file[0]
        ps_kwargs["h5_parent_group_name"] = self.cur_h5_file[1]

        ps_kwargs["optimize_function"] = self.optimize_function
        ps_kwargs["optimize_kwargs"] = self.combined_optimize_defaults

        ps_kwargs["reinitialize_function"] = self.initialize_function
        ps_kwargs["reinitialize_kwargs"] = self.combined_init_defaults
        ps_kwargs["reinitialize_before_sweep"] = self.reinitialize_before_sweep

        ps_kwargs["custom_do_param_sweep"] = self.custom_do_param_sweep
        ps_kwargs["custom_do_param_sweep_kwargs"] = self.custom_do_param_sweep_kwargs

        ps_kwargs["probe_function"] = self.probe_function

        # ps_kwargs["number_of_subprocesses"] = self.number_of_subprocesses
        # LEGACY will need to change when parallmanage radded to diff
        ps_kwargs[
            "differential_sweep_specs"
        ] = ParameterSweepReader()._dict_to_diff_spec(m, self.differential_sweep_specs)
        ps = DifferentialParameterSweep(**ps_kwargs)

        # these are added for custom function access
        ps.build_model = self.build_function
        ps.build_sweep_params = ParameterSweepReader()._dict_to_params
        ps.build_outputs = None
        ps.build_model_kwargs = self.combined_build_defaults
        ps.build_sweep_params_kwargs = {"input_dict": self.sweep_params}

        ps.parameter_sweep(
            m,
            ParameterSweepReader()._dict_to_params(m, input_dict=self.sweep_params),
            num_samples=self.num_samples,
            # build_model_kwargs=self.combined_build_defaults,
            # build_sweep_params_kwargs={"input_dict": self.sweep_params},
        )
