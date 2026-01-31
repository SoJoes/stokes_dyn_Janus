# REQUIRES PYTENTIAL

# based off of code in https://github.com/inducer/sumpy/blob/main/sumpy/kernel.py#L483-L524

import sumpy.kernel
from sumpy.kernel import ExpressionKernel, KernelArgument
from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    ClassVar,
    Generic,
    Literal,
    TypeVar,
    cast,
    overload,
)

import numpy as np
from typing_extensions import override

import loopy as lp
import pymbolic.primitives as prim
from pymbolic import Expression, var
from pymbolic.mapper import CSECachingMapperMixin, IdentityMapper
from pymbolic.primitives import make_sym_vector
from pytools import memoize_method

import sumpy.symbolic as symp
from sumpy.derivative_taker import (
    DerivativeCoeffDict,
    DifferentiatedExprDerivativeTaker,
    ExprDerivativeTaker,
    diff_derivative_coeff_dict,
)
from sumpy.symbolic import SpatialConstant, pymbolic_real_norm_2
from sumpy.derivative_taker import RadialDerivativeTaker


from collections.abc import Callable, Iterable, Sequence

import sympy as sp

from sumpy.assignment_collection import SymbolicAssignmentCollection
from sumpy.expansion.diff_op import LinearPDESystemOperator

@dataclass(frozen=True, repr=False)
class ScreenedLaplaceKernel(ExpressionKernel):
    """
    .. autoattribute:: lambda_name
    .. autoattribute:: allow_evanescent
    """

    mapper_method: ClassVar[str] = "map_expression_kernel" # TODO: CHOOSE THIS

    lambda_name: str
    """The argument name to use for the Screened Laplace parameter when generating
    functions to evaluate this kernel.
    """
    allow_evanescent: bool

    def __init__(self,
                 dim: int,
                 lambda_name: str = "lam",
                 allow_evanescent: bool = False) -> None:
        lam = SpatialConstant(lambda_name)

        # Guard against code using the old positional interface.
        assert isinstance(allow_evanescent, bool)

        if dim == 2:
            r = pymbolic_real_norm_2(make_sym_vector("d", dim))
            expr = var("hankel_1")(0, var("I")*lam*r)
            scaling = var("pi")/2*var("I")
        else:
            raise NotImplementedError(f"unsupported dimension: '{dim}'")

        super().__init__(dim, expression=expr, global_scaling_const=scaling)

        object.__setattr__(self, "lambda_name", lambda_name)
        object.__setattr__(self, "allow_evanescent", allow_evanescent)

    @property
    @override
    def is_complex_valued(self) -> bool:
        return True

    @override
    def __str__(self) -> str:
        return f"ScrnLplKnl{self.dim}D({self.lambda_name})"

    @override
    def prepare_loopy_kernel(self, loopy_knl: lp.TranslationUnit) -> lp.TranslationUnit:
        from sumpy.codegen import register_bessel_callables
        return register_bessel_callables(loopy_knl)

    @override
    def get_args(self) -> Sequence[KernelArgument]:
        k_dtype = np.complex128 if self.allow_evanescent else np.float64
        return [
                KernelArgument(
                    loopy_arg=lp.ValueArg(self.lambda_name, k_dtype),
                    )]

    @override
    def get_derivative_taker(
            self,
            dvec: symp.Matrix,
            rscale: symp.Expr,
            sac: SymbolicAssignmentCollection | None,
        ) -> ExprDerivativeTaker:
        from sumpy.derivative_taker import RadialDerivativeTaker
        return RadialDerivativeTaker(self.get_expression(dvec), dvec, rscale, sac)

    @override
    def get_pde_as_diff_op(self) -> LinearPDESystemOperator:
        from sumpy.expansion.diff_op import laplacian, make_identity_diff_op

        w = make_identity_diff_op(self.dim)
        lam = symp.Symbol(self.lambda_name)
        return laplacian(w) - lam**2 * w
