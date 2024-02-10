"""Stolen from https://github.com/ComputationalPhysiology/simcardems/blob/main/src/simcardems/utils.py"""

import dolfin
import ufl_legacy as ufl
import numpy as np


class Projector:
    def __init__(
        self,
        V: dolfin.FunctionSpace,
        solver_type: str = "lu",
        preconditioner_type: str = "default",
    ):
        """
        Projection class caching solver and matrix assembly

        Args:
            V (dolfin.FunctionSpace): Function-space to project in to
            solver_type (str, optional): Type of solver. Defaults to "lu".
            preconditioner_type (str, optional): Type of preconditioner. Defaults to "default".

        Raises:
            RuntimeError: _description_
        """
        u = dolfin.TrialFunction(V)
        self._v = dolfin.TestFunction(V)
        self._dx = dolfin.Measure("dx", domain=V.mesh())
        self._b = dolfin.Function(V)
        self._A = dolfin.assemble(ufl.inner(u, self._v) * self._dx)
        lu_methods = dolfin.lu_solver_methods().keys()
        krylov_methods = dolfin.krylov_solver_methods().keys()
        if solver_type == "lu" or solver_type in lu_methods:
            if preconditioner_type != "default":
                raise RuntimeError("LUSolver cannot be preconditioned")
            self.solver = dolfin.LUSolver(self._A, "default")
        elif solver_type in krylov_methods:
            self.solver = dolfin.PETScKrylovSolver(
                self._A,
                solver_type,
                preconditioner_type,
            )
        else:
            raise RuntimeError(
                f"Unknown solver type: {solver_type}, method has to be lu"
                + f", or {np.hstack(lu_methods, krylov_methods)}",
            )
        self.solver.set_operator(self._A)

    def project(self, u: dolfin.Function, f: ufl.core.expr.Expr) -> None:
        """
        Project `f` into `u`.

        Args:
            u (dolfin.Function): The function to project into
            f (ufl.core.expr.Expr): The ufl expression to project
        """
        dolfin.assemble(ufl.inner(f, self._v) * self._dx, tensor=self._b.vector())
        self.solver.solve(u.vector(), self._b.vector())

    def __call__(self, u: dolfin.Function, f: ufl.core.expr.Expr) -> None:
        self.project(u=u, f=f)


def von_mises(T: ufl.Coefficient) -> ufl.Coefficient:
    r"""Compute the von Mises stress tensor :math`\sigma_v`, with

    .. math::

        \sigma_v^2 = \frac{1}{2} \left(
            (\mathrm{T}_{11} - \mathrm{T}_{22})^2 +
            (\mathrm{T}_{22} - \mathrm{T}_{33})^2 +
            (\mathrm{T}_{33} - \mathrm{T}_{11})^2 +
        \right) - 3 \left(
            \mathrm{T}_{12} + \mathrm{T}_{23} + \mathrm{T}_{31}
        \right)

    Parameters
    ----------
    T : ufl.Coefficient
        Cauchy stress tensor

    Returns
    -------
    ufl.Coefficient
        The von Mises stress tensor
    """
    von_Mises_squared = 0.5 * (
        (T[0, 0] - T[1, 1]) ** 2 + (T[1, 1] - T[2, 2]) ** 2 + (T[2, 2] - T[0, 0]) ** 2
    ) + 3 * (T[0, 1] + T[1, 2] + T[2, 0])

    return ufl.sqrt(abs(von_Mises_squared))


def ca_transient(t, tstart=0.05, ca_ampl=0.3):
    tau1 = 0.05
    tau2 = 0.110

    ca_diast = 0.0

    beta = (tau1 / tau2) ** (-1 / (tau1 / tau2 - 1)) - (tau1 / tau2) ** (-1 / (1 - tau2 / tau1))
    ca = np.zeros_like(t)

    ca[t <= tstart] = ca_diast

    ca[t > tstart] = (ca_ampl - ca_diast) / beta * (
        np.exp(-(t[t > tstart] - tstart) / tau1) - np.exp(-(t[t > tstart] - tstart) / tau2)
    ) + ca_diast
    return ca


def try_except_runtimererror(func, values, args):
    try:
        values.append(func(*args))
    except RuntimeError:
        pass


def name2latex(name: str) -> str:
    name_lst = name.split("_")
    if len(name_lst) == 1:
        return f"${name}$"

    *sym, sup = name_lst
    if len(sym) == 1:
        sym0 = sym[0]
        if sym0 == "sigma":
            sym0 = "\\sigma"
        return f"${sym0}_{{{sup}}}$"
    else:
        if sym[0] == "sigma":
            sym0 = "\\mathrm{dev} \\sigma"
            return f"${sym0}_{{{sup}}}$"

        raise ValueError(sym)
