import subprocess as sp
from pathlib import Path


results_folder_tmp = (
    "/Users/finsberg/local/src/ciucci-2024/code/results/cylinder_{varying_str}/spring3000.0"
)
figure_folder_tmp = (
    "/Users/finsberg/local/src/ciucci-2024/figures/cylinder_{varying_str}/spring3000.0"
)


def create_uniform_gamma_figures():
    print("Uniform gamma")

    figure_folder = figure_folder_tmp.format(varying_str="")
    results_folder = results_folder_tmp.format(varying_str="")

    Path(figure_folder).mkdir(exist_ok=True, parents=True)

    for min_value, max_value, name, label, blocknr in [
        (-1.0, 1.0, "sigma_xx", "$\\sigma_{xx}$", 1),
        (-1.0, 1.0, "sigma_r", "$\\sigma_{r}$", 2),
        (-1.0, 1.0, "sigma_c", "$\\sigma_{c}$", 3),
        (-1.0, 1.0, "sigma_dev_xx", "dev $\\sigma_{xx}$", 4),
        (-1.0, 1.0, "sigma_dev_r", "dev $\\sigma_{r}$", 5),
        (-1.0, 1.0, "sigma_dev_c", "dev $\\sigma_{c}$", 6),
        (-0.3, 0.0, "E_xx", "$E_{xx}$", 7),
        (0.0, 0.2, "E_r", "$E_{r}$", 8),
        (0.0, 0.2, "E_c", "$E_{c}$", 9),
        (-1.0, 1.0, "hydrostatic", "$p$", 10),
    ]:
        print(name)
        sp.run(
            [
                "/Applications/ParaView-5.11.0.app/Contents/bin/pvpython",
                "create_paraview_figure.py",
                "--filename",
                f"{results_folder}/results.xdmf",
                "--figure_folder",
                figure_folder,
                "--min_value",
                str(min_value),
                "--max_value",
                str(max_value),
                "--name",
                name,
                "--label",
                label,
                "--blocknr",
                str(blocknr),
            ]
        )


def create_varing_gamma_figures():
    print("Varying gamma")
    figure_folder = figure_folder_tmp.format(varying_str="_varying")
    results_folder = results_folder_tmp.format(varying_str="_varying")

    Path(figure_folder).mkdir(exist_ok=True, parents=True)

    for min_value, max_value, name, label, blocknr in [
        (-1.0, 2.0, "sigma_xx", "$\\sigma_{xx}$", 1),
        (-1.0, 1.0, "sigma_r", "$\\sigma_{r}$", 2),
        (-1.0, 1.0, "sigma_c", "$\\sigma_{c}$", 3),
        (-1.0, 2.0, "sigma_dev_xx", "dev $\\sigma_{xx}$", 4),
        (-2.0, 1.0, "sigma_dev_r", "dev $\\sigma_{r}$", 5),
        (-2.0, 1.0, "sigma_dev_c", "dev $\\sigma_{c}$", 6),
        (-0.3, 0.0, "E_xx", "$E_{xx}$", 7),
        (0.0, 0.2, "E_r", "$E_{r}$", 8),
        (0.0, 0.2, "E_c", "$E_{c}$", 9),
        (-2.0, 1.0, "hydrostatic", "$p$", 10),
    ]:
        print(name)

        sp.run(
            [
                "/Applications/ParaView-5.11.0.app/Contents/bin/pvpython",
                "create_paraview_figure.py",
                "--filename",
                f"{results_folder}/results.xdmf",
                "--figure_folder",
                figure_folder,
                "--min_value",
                str(min_value),
                "--max_value",
                str(max_value),
                "--name",
                name,
                "--label",
                label,
                "--blocknr",
                str(blocknr),
            ]
        )


def main():
    # create_uniform_gamma_figures()
    create_varing_gamma_figures()


if __name__ == "__main__":
    main()
