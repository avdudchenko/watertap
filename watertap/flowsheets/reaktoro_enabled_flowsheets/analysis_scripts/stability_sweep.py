import platform
import os
from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheets.softening_acid_ro as sar


def main():
    cwd = get_working_dir()
    loopTool(
        cwd + "/stability_sweep_config.yaml",
        build_function=sar.build_model,
        initialize_function=sar.initialize,
        optimize_function=sar.solve_model,
        save_name="stability_sweep",
        probe_function=sar.test_func,
        saving_dir=cwd,
        number_of_subprocesses=1,
        num_loop_workers=3,
    )


if __name__ == "__main__":
    main()
