from typing import Literal
from dataclasses import dataclass
import shutil
import dolfin
from pathlib import Path

import pulse

import numpy as np
import ufl_legacy as ufl
import utils

from geometry import get_lv_geometry


dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True


class SmoothLV(dolfin.UserExpression):
    def __init__(self, f):
        self.f = f
        super().__init__()

    def eval(self, value, x):
        dx = 0.05
        values = [self.f(x[0], x[1], x[2])]

        for i in [-2, -1, 1, 2]:
            utils.try_except_runtimererror(self.f, values, (x[0] + i * dx, x[1], x[2]))
            utils.try_except_runtimererror(self.f, values, (x[0], x[1] + i * dx, x[2]))
            utils.try_except_runtimererror(self.f, values, (x[0], x[1], x[2] + i * dx))

        value[0] = np.mean(values)

    def value_shape(self):
        return ()


@dataclass
class DataCollector:
    problem: pulse.MechanicsProblem
    folder: Path

    @property
    def geo(self):
        return self.problem.geometry

    @property
    def material(self):
        return self.problem.material

    def __post_init__(self):
        shutil.rmtree(self.folder, ignore_errors=True)
        self.folder.mkdir(exist_ok=True, parents=True)
        self.path_reference = self.folder / "results_reference.xdmf"
        self.path_current = self.folder / "results_current.xdmf"
        # self.path_reference_smooth = self.folder / "results_reference_smooth.xdmf"
        # self.path_current_smooth = self.folder / "results_current_smooth.xdmf"

        self.functions = {
            "reference": {},
            "current": {},
            # "reference_smooth": {},
            # "current_smooth": {},
        }

        self.reference_mesh = self.geo.mesh
        self.current_mesh = dolfin.Mesh(self.geo.mesh)

        self.V_int = dolfin.VectorFunctionSpace(self.current_mesh, "CG", 1)
        self.u_int = dolfin.Function(self.V_int)

        for mesh, label in [(self.reference_mesh, "reference"), (self.current_mesh, "current")]:
            self.functions[label]["u"] = dolfin.Function(
                dolfin.VectorFunctionSpace(mesh, "CG", 2), name="u"
            )
            self.functions[label]["p"] = dolfin.Function(
                dolfin.FunctionSpace(mesh, "CG", 1), name="u"
            )

            V_DG2 = dolfin.FunctionSpace(mesh, "DG", 0)
            self.functions[label]["sigma_ff"] = dolfin.Function(V_DG2)
            self.functions[label]["sigma_ss"] = dolfin.Function(V_DG2)
            self.functions[label]["sigma_nn"] = dolfin.Function(V_DG2)
            self.functions[label]["sigma_dev_ff"] = dolfin.Function(V_DG2)
            self.functions[label]["sigma_dev_ss"] = dolfin.Function(V_DG2)
            self.functions[label]["sigma_dev_nn"] = dolfin.Function(V_DG2)
            self.functions[label]["E_ff"] = dolfin.Function(V_DG2)
            self.functions[label]["E_ss"] = dolfin.Function(V_DG2)
            self.functions[label]["E_nn"] = dolfin.Function(V_DG2)
            self.functions[label]["von_Mises"] = dolfin.Function(V_DG2)

            # V_smooth = dolfin.FunctionSpace(mesh, "DG", 0)
            # self.functions[f"{label}_smooth"]["sigma_ff"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["sigma_ss"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["sigma_nn"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["sigma_dev_ff"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["sigma_dev_ss"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["sigma_dev_nn"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["E_ff"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["E_ss"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["E_nn"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["von_Mises"] = dolfin.Function(V_smooth)
            # self.functions[f"{label}_smooth"]["p"] = dolfin.Function(V_smooth)

        self.proj = utils.Projector(self.functions["reference"]["sigma_ff"].function_space())

    def save(self, ti: float):
        u, p = self.problem.state.split(deepcopy=True)
        F = pulse.kinematics.DeformationGradient(u)
        Fe = self.material.Fe(F)
        sigma = self.material.CauchyStress(F, p)
        sigma_dev = sigma - (1 / 3) * ufl.tr(sigma) * ufl.Identity(3)
        E = pulse.kinematics.GreenLagrangeStrain(Fe)
        f = F * self.geo.f0
        s = F * self.geo.s0
        n = F * self.geo.n0

        # First compute values on the reference mesh
        self.functions["reference"]["u"].assign(u)
        self.functions["reference"]["p"].assign(p)
        self.proj.project(self.functions["reference"]["sigma_ff"], dolfin.inner(f, sigma * f))
        self.proj.project(self.functions["reference"]["sigma_ss"], dolfin.inner(s, sigma * s))
        self.proj.project(self.functions["reference"]["sigma_nn"], dolfin.inner(n, sigma * n))
        self.proj.project(
            self.functions["reference"]["sigma_dev_ff"], dolfin.inner(f, sigma_dev * f)
        )
        self.proj.project(
            self.functions["reference"]["sigma_dev_ss"], dolfin.inner(s, sigma_dev * s)
        )
        self.proj.project(
            self.functions["reference"]["sigma_dev_nn"], dolfin.inner(n, sigma_dev * n)
        )
        self.proj.project(
            self.functions["reference"]["E_ff"], dolfin.inner(self.geo.f0, E * self.geo.f0)
        )
        self.proj.project(
            self.functions["reference"]["E_ss"], dolfin.inner(self.geo.s0, E * self.geo.s0)
        )
        self.proj.project(
            self.functions["reference"]["E_nn"], dolfin.inner(self.geo.n0, E * self.geo.n0)
        )
        self.proj.project(self.functions["reference"]["von_Mises"], utils.von_mises(sigma))

        # Then transfer to the current mesh
        for name, f in self.functions["current"].items():
            f.vector()[:] = self.functions["reference"][name].vector()[:]

        # Move mesh
        u_int = dolfin.interpolate(self.functions["current"]["u"], self.V_int)
        # Move mesh back
        self.u_int.vector()[:] = -self.u_int.vector()[:]
        dolfin.ALE.move(self.current_mesh, self.u_int)
        # Move mesh to new position
        dolfin.ALE.move(self.current_mesh, u_int)
        self.u_int.vector()[:] = u_int.vector()[:]

        # Save to file
        for label, path in [
            ("reference", self.path_reference),
            ("current", self.path_current),
        ]:
            with dolfin.XDMFFile(path.as_posix()) as xdmf:
                if label in ["reference", "current"]:
                    xdmf.write_checkpoint(
                        self.functions[label]["u"],
                        function_name="u",
                        time_step=ti,
                        encoding=dolfin.XDMFFile.Encoding.HDF5,
                        append=True,
                    )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_ff"],
                    function_name="sigma_ff",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_ss"],
                    function_name="sigma_ss",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_nn"],
                    function_name="sigma_nn",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_dev_ff"],
                    function_name="sigma_dev_ff",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_dev_ss"],
                    function_name="sigma_dev_ss",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["sigma_dev_nn"],
                    function_name="sigma_dev_nn",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["E_ff"],
                    function_name="E_ff",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["E_ss"],
                    function_name="E_ss",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["E_nn"],
                    function_name="E_nn",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["von_Mises"],
                    function_name="von_Mises",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )
                xdmf.write_checkpoint(
                    self.functions[label]["p"],
                    function_name="p",
                    time_step=ti,
                    encoding=dolfin.XDMFFile.Encoding.HDF5,
                    append=True,
                )


def main(
    output_folder,
    mesh_folder: Path = Path("meshes/lv"),
    case: Literal["native", "transplanted"] = "native",
):
    geo = get_lv_geometry(mesh_folder=mesh_folder)

    microstructure = pulse.Microstructure(f0=geo.f0, s0=geo.s0, n0=geo.n0)

    geometry = pulse.HeartGeometry(
        mesh=geo.mesh,
        markers=geo.markers,
        marker_functions=pulse.MarkerFunctions(ffun=geo.ffun),
        microstructure=microstructure,
    )

    m2mm = 1000.0

    matparams = {
        "a": 2.280,
        "b": 9.726,
        "a_f": 1.685,
        "b_f": 15.779,
        "a_s": 0.0,
        "b_s": 0.0,
        "a_fs": 0.0,
        "b_fs": 0.0,
    }
    gamma = dolfin.Constant(0.0)
    activation = gamma

    material = pulse.HolzapfelOgden(
        active_model="active_strain",
        activation=activation,
        parameters=matparams,
        f0=geometry.f0,
        s0=geometry.s0,
        n0=geometry.n0,
    )

    # Pericardium type Robin BC
    spring = dolfin.Constant(500 / m2mm)  # kPa/mm
    robin_bc = [
        pulse.RobinBC(value=dolfin.Constant(spring), marker=geo.markers["EPI"][0]),
        # pulse.RobinBC(value=dolfin.Constant(spring), marker=geo.markers["BASE"][0]),
    ]

    # LV Pressure
    lvp = dolfin.Constant(0.0)
    lv_marker = geometry.markers["ENDO"][0]
    lv_pressure = pulse.NeumannBC(traction=lvp, marker=lv_marker, name="lv")
    neumann_bc = [lv_pressure]

    # Fix the basal plane in the longitudinal direction
    def fix_basal_plane(W):
        V = W if W.sub(0).num_sub_spaces() == 0 else W.sub(0)
        bc = dolfin.DirichletBC(
            V.sub(0),
            dolfin.Constant(0.0),
            geometry.ffun,
            geometry.markers["BASE"][0],
        )
        return bc

    dirichlet_bc = [fix_basal_plane]

    bcs = pulse.BoundaryConditions(dirichlet=dirichlet_bc, neumann=neumann_bc, robin=robin_bc)

    problem = pulse.MechanicsProblem(geometry, material, bcs)

    data_collector = DataCollector(problem=problem, folder=output_folder)

    # Pressure alues are taken from https://journals.physiology.org/doi/full/10.1152/ajpheart.00218.2013
    # Volumes are given in spreadsheet from data
    if case == "native":
        EDP = 1.3
        ESP = 13.5
        target_EDV = 115.8
        target_ESV = 52.2
        gamma_ES = 0.303
        # gamma_ES = 124.0

    elif case == "transplanted":
        EDP = 1.0
        ESP = 8.0
        target_EDV = 37.4
        target_ESV = 34.1
        # gamma_ES = 0.165
        gamma_ES = 0.168
        # gamma_ES = 39.5

    gammas = [0.0, 0.0, gamma_ES]
    pressures = [0.0, EDP, ESP]
    volumes = [geometry.cavity_volume(u=problem.state.split()[0])]

    np.save(output_folder / "gammas.npy", gammas)
    np.save(output_folder / "pressures.npy", pressures)

    data_collector.save(0.0)

    # ED
    pulse.iterate.iterate(
        problem,
        lvp,
        EDP,
        initial_number_of_steps=20,
        continuation=False,
    )
    EDV = geometry.cavity_volume(u=problem.state.split()[0])
    volumes.append(EDV)
    data_collector.save(1.0)
    print(f"EDP: {EDP} kPa, EDV: {EDV} uL, target EDV: {target_EDV}, gamma: {float(gamma)}")

    # ES
    pulse.iterate.iterate(
        problem,
        (lvp, gamma),
        (ESP, gamma_ES),
        initial_number_of_steps=1000,
        continuation=True,
    )
    ESV = geometry.cavity_volume(u=problem.state.split()[0])
    volumes.append(ESV)
    print(f"ESP: {ESP} kPa, ESV: {ESV} uL, target ESV: {target_ESV}, gamma: {float(gamma)}")
    data_collector.save(2.0)

    np.save(output_folder / "volumes.npy", volumes)
