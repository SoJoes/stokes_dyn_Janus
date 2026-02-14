import numpy as np

from meshmode.array_context import PyOpenCLArrayContext
from meshmode.discretization import Discretization
from meshmode.discretization.poly_element import (
    InterpolatoryQuadratureSimplexGroupFactory,
)
from meshmode.dof_array import DOFArray

from pytential import bind, sym
from pytential.target import PointsTarget

# my file
from pytential_handling.my_laplace_kernel import ScreenedLaplaceKernel
from settings import actx
import sys
from os import devnull


# Disable
def blockPrint():
    sys.stdout = open(devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

# {{{ set some constants for use below

nelements = 20
bdry_quad_order = 4 # order of quadrature on the boundary
mesh_order = bdry_quad_order
qbx_order = bdry_quad_order
bdry_ovsmp_quad_order = 4*bdry_quad_order # boundary ? quadrature order
fmm_order = 10
k = 0.1

# }}}

from meshmode.mesh.generation import ellipse, make_curve_mesh
from functools import partial
from meshmode.mesh.processing import affine_map, merge_disjoint_meshes

class Janus_particle_array:
  def __init__(self, positions, facings, base_mesh, mesh_scale=1):
    self.pos_array = positions
    self.facing_array = facings
    self.base_mesh = base_mesh
    self.meshes = [affine_map(
                    self.base_mesh,
                    A=np.diag([mesh_scale, mesh_scale]),
                    b=position) for position in positions]

def amphilics(visualize=False, particle_pos=None, particle_facing=None):
    base_mesh = make_curve_mesh(
                partial(ellipse, 1),
                np.linspace(0, 1, nelements+1),
                mesh_order)

    Januses = Janus_particle_array(
        positions = particle_pos,
        facings = particle_facing,
        base_mesh = base_mesh
    )

    meshes = Januses.meshes

    mesh = merge_disjoint_meshes(meshes, single_group=False) # so that we have separate particles

    if visualize:
        from meshmode.mesh.visualization import draw_curve
        draw_curve(mesh)
        import matplotlib.pyplot as plt
        plt.show()

    # discretise boundary before qbx or gmres
    pre_density_discr = Discretization(
            actx, mesh,
            InterpolatoryQuadratureSimplexGroupFactory(bdry_quad_order))

    from pytential.qbx import QBXLayerPotentialSource, QBXTargetAssociationFailedError
    # setup for qbx
    qbx = QBXLayerPotentialSource(
            pre_density_discr,
            fine_order=bdry_ovsmp_quad_order,
            qbx_order=qbx_order,
            fmm_order=fmm_order
            )

    from sumpy.visualization import FieldPlotter
    fplot = FieldPlotter(np.zeros(2), extent=5, npoints=500)
    targets = actx.from_numpy(fplot.points)

    from pytential import GeometryCollection
    places = GeometryCollection({
        "qbx": qbx,
        "qbx_high_target_assoc_tol": qbx.copy(target_association_tolerance=0.05),
        "targets": PointsTarget(targets)
        }, auto_where="qbx") # places in which to eval qbx

    # discretised boundary for density calculations
    density_discr = places.get_discretization("qbx")

    # {{{ describe bvp

    from sumpy.kernel import HelmholtzKernel, LaplaceKernel
    kernel = ScreenedLaplaceKernel(2)

    sigma_sym = sym.var("sigma")
    sqrt_w = sym.sqrt_jac_q_weight(2)
    inv_sqrt_w_sigma = sym.cse(sigma_sym/sqrt_w)

    loc_sign = -1 #exterior condition DO NOT CHANGE

    k_sym = sym.var("k")
    bdry_op_sym = (-loc_sign*0.5*sigma_sym
            + sqrt_w*(
                sym.S(kernel, inv_sqrt_w_sigma, lam=k_sym,
                    qbx_forced_limit=+1)
                + sym.D(kernel, inv_sqrt_w_sigma, lam=k_sym,
                    qbx_forced_limit="avg")
                ))

    # }}}


    bound_op = bind(places, bdry_op_sym)

    # {{{ fix rhs and solve

    nodes = actx.thaw(density_discr.nodes())

    def amphilic(nodes, janus_array):
      x, y = nodes

      bc_data = []

      for igrp, facing in enumerate(janus_array.facing_array):
          cos_f = actx.np.cos(facing)
          sin_f = actx.np.sin(facing)

          # change center of rotation to center of particle
          xg = x[igrp] - janus_array.pos_array[igrp,0]
          yg = y[igrp] - janus_array.pos_array[igrp,1]

          rot_x =  cos_f * xg + sin_f * yg
          rot_y = -sin_f * xg + cos_f * yg

          theta = actx.np.arctan2(rot_y, rot_x)
          bc_g = (actx.np.cos(theta) + 1) / 2

          bc_data.append(bc_g)

      return DOFArray(actx, tuple(bc_data))

    bc = amphilic(nodes, Januses)

    bvp_rhs = bind(places, sqrt_w*sym.var("bc"))(actx, bc=bc)


    bvp_rhs = DOFArray(actx, data=tuple(array for array in bvp_rhs))

    from pytential.linalg.gmres import gmres

    blockPrint()
    gmres_result = gmres(
            bound_op.scipy_op(actx, sigma_sym.name, dtype=np.complex128, k=k),
            bvp_rhs, tol=1e-8, progress=True,
            stall_iterations=0,
            hard_failure=True) # figure out gmres
    enablePrint()

    print("gmres succeeded")

    # }}}

    # {{{ postprocess/visualize

    repr_kwargs = {
            "source": "qbx_high_target_assoc_tol",
            "target": "targets",
            "qbx_forced_limit": None}
    representation_sym = (
            sym.S(kernel, inv_sqrt_w_sigma,lam=k_sym, **repr_kwargs)
            + sym.D(kernel, inv_sqrt_w_sigma, lam=k_sym, **repr_kwargs))
    ones_density = density_discr.zeros(actx)
    for elem in ones_density:
        elem.fill(1)

    indicator = actx.to_numpy(
            bind(places, sym.D(LaplaceKernel(2), sigma_sym, **repr_kwargs))(
                actx, sigma=ones_density))

    try:
        fld_in_vol = actx.to_numpy(
                bind(places, representation_sym)(
                    actx, sigma=gmres_result.solution, k=k)).astype(np.float64)
    except QBXTargetAssociationFailedError as e:
        fplot.write_vtk_file("BIE-Janus-failed-targets.vts", [
            ("failed", actx.to_numpy(e.failed_target_flags))
            ])
        raise

    # calculate force on particle

    # find grad of potential
    from pytential.symbolic.primitives import grad

    representation_sym_grad = -grad(ambient_dim=2, operand = representation_sym)  # gives wrong result if not negative
    nabla_pot = bind(places, representation_sym_grad)(
      actx, sigma=gmres_result.solution, k=k)

    # nabla_pot is a vector field (2 components), so write components separately to VTK
    nabla_pot_x = actx.to_numpy(nabla_pot[0])
    nabla_pot_y = actx.to_numpy(nabla_pot[1])

    # calculate hydrophobic stress
    def hydrophobic_stress_T(u_sym, grad_u_sym, gamma=1, rho=1):
      # grad_u_sym is expected to be a symbolic vector (e.g., a tuple of expressions)
      grad_x_sym = grad_u_sym[0]
      grad_y_sym = grad_u_sym[1]

      # Magnitude squared of gradient
      grad_mag_sq = grad_x_sym**2 + grad_y_sym**2

      # Scalar part of the first two terms in the definition of T_ij
      # (u^2/rho) * delta_ij + (1/2) * |grad u|^2 * delta_ij
      scalar_diagonal_term = (u_sym**2 / rho) + rho * (1/2) * grad_mag_sq

      # Factor for the outer product term: -2 * rho * (grad_i u) * (grad_j u)
      outer_product_factor = -2 * rho

      T_xx_sym = gamma * (scalar_diagonal_term + outer_product_factor * grad_x_sym * grad_x_sym)
      T_xy_sym = gamma * (outer_product_factor * grad_x_sym * grad_y_sym)
      T_yx_sym = T_xy_sym
      T_yy_sym = gamma * (scalar_diagonal_term + outer_product_factor * grad_y_sym * grad_y_sym)

      # Return components of the stress tensor as a tuple
      return (T_xx_sym, T_xy_sym, T_yx_sym, T_yy_sym)

    T_sym_components = hydrophobic_stress_T(representation_sym, representation_sym_grad)

    # Get normal vectors for the density discretization
    mv_normal = bind(density_discr, sym.normal(2))(actx)
    normal = mv_normal.as_vector(object)

    nvec_sym = sym.make_sym_vector("normal", 2)

    # Define a symbolic representation for the potential *on the boundary itself*
    repr_kwargs_boundary = {
            "source": "qbx",  # Source is 'qbx' (pre_density_discr)
            "target": "qbx",  # Target is 'qbx' (the boundary itself)
            "qbx_forced_limit": +1 # Or appropriate limit for boundary evaluation
    }
    representation_sym_boundary = (
            sym.S(kernel, inv_sqrt_w_sigma, lam=k_sym, **repr_kwargs_boundary)
            + sym.D(kernel, inv_sqrt_w_sigma, lam=k_sym, **repr_kwargs_boundary)
    )

    # find grad of potential on the boundary
    representation_sym_grad_boundary = grad(ambient_dim=2, operand = representation_sym_boundary)

    # calculate hydrophobic stress tensor on the boundary
    T_sym_components_boundary = hydrophobic_stress_T(representation_sym_boundary, representation_sym_grad_boundary)

    # Define force integrands
    force_integrand_x_sym = T_sym_components_boundary[0] * nvec_sym[0] + T_sym_components_boundary[1] * nvec_sym[1]
    force_integrand_y_sym = T_sym_components_boundary[2] * nvec_sym[0] + T_sym_components_boundary[3] * nvec_sym[1]

    # TORQUE
    mv_pos = bind(density_discr, sym.nodes(2))(actx)
    pos = mv_pos.as_vector(object)

    # Use separate symbolic variables for position components for evaluation compatibility
    r_pos_sym = sym.make_sym_vector("r_pos", 2)

    # formula derived from definitions by me <3
    torque_integrand_sym = r_pos_sym[0] * force_integrand_y_sym - r_pos_sym[1] * force_integrand_x_sym

    # Use sym.integral to sum the integrands over the boundary
    from pytential.symbolic.primitives import integral as sym_integral
    from pytential.symbolic.primitives import area_element, QWeight

    force_density_x = bind(places, force_integrand_x_sym)(actx,sigma=gmres_result.solution,k=k,normal=normal)
    force_density_y = bind(places,force_integrand_y_sym)(actx,sigma=gmres_result.solution,k=k,normal=normal)
    torque_density = bind(places, torque_integrand_sym)(actx, sigma=gmres_result.solution, k=k, normal=normal, r_pos=pos)

    dS = area_element(1, 1, None) * QWeight(None)

    integral_weights = bind(density_discr, dS)(actx)

    n_particles = len(force_density_x)

    forces_x = np.ones(n_particles)
    forces_y = np.ones(n_particles)
    torques = np.ones(n_particles)

    for igrp in range(n_particles):
        # manual node.sum
        fx = actx.to_numpy(
            actx.np.sum(force_density_x[igrp] * integral_weights[igrp])
          )
        fy = actx.to_numpy(
            actx.np.sum(force_density_y[igrp] * integral_weights[igrp])
          )
        t = actx.to_numpy(
            actx.np.sum(torque_density[igrp] * integral_weights[igrp])
            )

        forces_x[igrp] = fx
        forces_y[igrp] = fy

        # torque needs to be centred around particle
        torques[igrp] = t - (Januses.pos_array[igrp][0] * fy - Januses.pos_array[igrp][1] * fx)

    return(forces_x, forces_y, torques)
