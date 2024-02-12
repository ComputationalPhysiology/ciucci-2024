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


def main(output_folder, mesh_folder: Path = Path("meshes/lv")):
    overwrite = True
    Path(output_folder).mkdir(exist_ok=True, parents=True)
    output = Path(output_folder) / "results.xdmf"
    if output.is_file() and not overwrite:
        print(f"Output {output} already exists")
        return

    geo = get_lv_geometry(mesh_folder=mesh_folder)

    V_DG2 = dolfin.FunctionSpace(geo.mesh, "DG", 1)
    proj = utils.Projector(V_DG2)
    sigma_ff = dolfin.Function(V_DG2)
    sigma_ss = dolfin.Function(V_DG2)
    sigma_nn = dolfin.Function(V_DG2)

    sigma_dev_ff = dolfin.Function(V_DG2)
    sigma_dev_ss = dolfin.Function(V_DG2)
    sigma_dev_nn = dolfin.Function(V_DG2)

    E_ff = dolfin.Function(V_DG2)
    E_ss = dolfin.Function(V_DG2)
    E_nn = dolfin.Function(V_DG2)

    von_Mises = dolfin.Function(V_DG2)

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
    gammas = [0.0, 0.2, 0.2]
    pressures = [0.0, 0.0, 15.0]

    np.save(output_folder / "gammas.npy", gammas)
    np.save(output_folder / "pressures.npy", pressures)

    for ti, (lvp_, g) in enumerate(zip(pressures, gammas)):
        print(lvp_, g)
        pulse.iterate.iterate(
            problem,
            (lvp, gamma),
            (lvp_, g),
            initial_number_of_steps=20,
            continuation=False,
        )
        u, p = problem.state.split()

        F = pulse.kinematics.DeformationGradient(u)
        sigma = material.CauchyStress(F, p)
        sigma_dev = sigma - (1 / 3) * ufl.tr(sigma) * ufl.Identity(3)
        E = pulse.kinematics.GreenLagrangeStrain(F)
        f = F * geo.f0
        s = F * geo.s0
        n = F * geo.n0

        proj.project(von_Mises, utils.von_mises(sigma))

        proj.project(sigma_ff, dolfin.inner(f, sigma * f))
        proj.project(sigma_ss, dolfin.inner(s, sigma * s))
        proj.project(sigma_nn, dolfin.inner(n, sigma * n))

        proj.project(sigma_dev_ff, dolfin.inner(f, sigma_dev * f))
        proj.project(sigma_dev_ss, dolfin.inner(s, sigma_dev * s))
        proj.project(sigma_dev_nn, dolfin.inner(n, sigma_dev * n))

        proj.project(E_ff, dolfin.inner(geo.f0, E * geo.f0))
        proj.project(E_ss, dolfin.inner(geo.s0, E * geo.s0))
        proj.project(E_nn, dolfin.inner(geo.n0, E * geo.n0))
        with dolfin.XDMFFile(output.as_posix()) as f:
            f.parameters["functions_share_mesh"] = True
            f.parameters["rewrite_function_mesh"] = False
            f.write_checkpoint(
                u,
                function_name="u",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_ff,
                function_name="sigma_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_ss,
                function_name="sigma_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_nn,
                function_name="sigma_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_ff,
                function_name="sigma_dev_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_ss,
                function_name="sigma_dev_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_nn,
                function_name="sigma_dev_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_ff,
                function_name="E_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_ss,
                function_name="E_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_nn,
                function_name="E_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                von_Mises,
                function_name="von_Mises",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                p,
                function_name="p",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
