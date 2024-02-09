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
