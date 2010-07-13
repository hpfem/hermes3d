.. _example-elasto-statics:

Elasto Statics
==============

**Git reference:** Examples `elasto statics <http://git.hpfem.org/hermes3d.git/tree/HEAD:/examples/elastostatics>`_.

This example describes the implementation of a linear elastic problem inside L-shape domain. Strictly speaking, 
elasticity problem cannot be solved exactly in general: they can be solved exactly in very particular cases 
(e.g., special geometry); The described problem in the example will be solved numerically. 

Elastostatics studys linear elastic problem under the conditions of equilibrium, where all forces on the elastic 
body sum to zero, and  displacements are not a function of time. 

.. index::
   single: mesh; fixed
   single: problem; elliptic, linear, symmetric

The governing equation under Cartesian coordinate:

.. math::
   :nowrap:
   :label: elasto-statics

   \begin{eqnarray*}
   \sigma_{ji,j} + F_i & = & 0 \hbox{ in }\Omega \\ \nonumber
   \epsilon_{ij}       & = & \frac{1}{2}(u_{j,i} + u_{i,j})   \\
   \sigma_{i,j}        & = & C_{ijkl} \, \epsilon_{kl}
   \end{eqnarray*}

where, the subscript $\cdot_{,j}$ indicate $\partial{\cdot}/\partial x_j$; $\sigma_{ji,j}$ is the 
stress tensor; $\epsilon_{ij}$ is the strain(deformation); $u_i$ is the displacement;
$C_{ijkl}$ is the forth order stiffness tensor; Please not that, by Einstein summation convention, 
the $3^{rd}$ equation of Equ :eq:`elasto-statics` represent the following: 

.. math::
   :nowrap:
   :label: elasto-sum

   \begin{eqnarray*}
   C_{ijkl} \, \epsilon_{kl} & = & \sum_{k,l=1}^3 C_{ijkl} \, \epsilon_{kl}
   \end{eqnarray*}

where $i, j, k, l$ all take on the values 1, 2, or 3. 

.. image:: elasto-statics-domain.png

Domain of interest is isotropic homogeneous media of L-shaped bean (see above), equipped with 
the zero Dirichlet boundary conditions: $u_1 = u_2 = u_3 = 0$ on the all 5 surfaces (${\Gamma}_u$) 
except the left vertical one (${\Gamma}_F$), where the external force $F$ is applied.  


.. code-block:: c++
::

        // Boundary condition types.
        BCType bc_types_x(int marker)
        {
          return BC_NATURAL;
        }

        BCType bc_types_y(int marker)
        {
          return BC_NATURAL;
        }

        BCType bc_types_z(int marker)
        {
          return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL;
        }


The stiffness tensor $C_{ijkl}$ is constant and symmetric due to homegeneity,  i.e., it has 
no preferred direction. Therefore, the stress tensor could be written as (derivation omitted):

.. math::
   :nowrap:
   :label: elasto-stress

   \begin{eqnarray*}
   \sigma_{ij} & = & \lambda \delta_{ij} \epsilon_{kk} + 2\mu\epsilon_{ij} \\ \nonumber
   \lambda     & = & \frac{E\nu}{(1+\nu)(1-2\nu)}                          \\
   \mu         & = & \frac{E}{2(1+\nu)} 
   \end{eqnarray*}

where $\lambda$ is the first lame parameter, $\mu$ is the second lame parameter or shear modulus, 
$E$ is the Young's modulus, $\nu$ is the Poisson's ratio. In our example, $E = 200 \times 10^9$ Gpa, 
$\nu = 0.3$. 

Substituting Equ :eq:`elasto-stress` back into Equ :eq:`elasto-statics` yield:
 
.. math::
   :nowrap:
   :label: elasto-navier

   \begin{eqnarray*}
   \mu u_{i,jj}  + (\mu + \lambda)u_{j,ij} + F_i & = & 0              \\ \nonumber
   \hbox{ or }           & \, & \\                                      
   \mu \Delta{u} + (\mu + \lambda) \mathsf{grad} \, \mathsf{div} u  + F & = & 0
   \end{eqnarray*}

The corresponding weak formulations are as following, each equation for displacement in one direction:

.. math::
   :nowrap:
   :label: elasto-statics-form

   \begin{eqnarray*}
   \int_{\Omega} (\lambda + 2\mu) u_{i} \, v_{i} + \mu u_{j} \, v_{j} + \mu u_{k} \, v_{k} \quad 
   +\quad \int_{\Omega} \lambda u_{i} \,  v_{j} + \mu u_{j} \, v_{i} \quad
   +\quad \int_{\Omega} \lambda u_{i} \,  v_{k} + \mu u_{k} \, v_{i}
     &  = & 0 \\ \nonumber
   \int_{\Omega} \mu u_{i} \, v_{i} + (\lambda + 2\mu) u_{j} \, v_{j} + \mu u_{k} \, v_{k} \quad
   +\quad \int_{\Omega} \lambda u_{j} \,  v_{k} + \mu u_{k} \, v_{j}
     &  = & 0 \\
   \int_{\Omega} \mu u_{i} \, v_{i} + \mu u_{j} \, v_{j} + (\lambda + 2\mu) u_{k} \, v_{k} 
     &  = & \int_{\Gamma_F} F_i v. \nonumber
   \end{eqnarray*}

Code for the weak forms:

.. code-block:: c++
::

        template<typename real, typename scalar>
        scalar bilinear_form_0_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dx_b_dy_c_dz<real, scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
        template<typename real, typename scalar>
        scalar bilinear_form_0_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dx_b_dy_c_dz<real, scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
        }

        template<typename real, typename scalar>
        scalar bilinear_form_0_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dudx_dvdy_b_dudy_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
        }

        template<typename real, typename scalar>
        scalar bilinear_form_0_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dudx_dvdz_b_dudz_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
        }

        template<typename real, typename scalar>
        scalar surf_linear_form_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return 0.0;
        }

        template<typename real, typename scalar>
        scalar bilinear_form_1_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dx_b_dy_c_dz<real, scalar>(mu, lambda + 2*mu, mu, n, wt, u, v, e);
        }

        template<typename real, typename scalar>
        scalar bilinear_form_1_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dudy_dvdz_b_dudz_dvdy<real, scalar>(lambda, mu, n, wt, v, u, e);
        }

        template<typename real, typename scalar>
        scalar surf_linear_form_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return 0.0;
        }

        template<typename real, typename scalar>
        scalar bilinear_form_2_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          return int_a_dx_b_dy_c_dz<real, scalar>(mu, mu, lambda + 2*mu, n, wt, u, v, e);
        }

        template<typename real, typename scalar>
        scalar surf_linear_form_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
        {
          scalar res = 0.0;
          for (int i = 0; i < n; i++)
            res += wt[i] * (f * v->fn[i]);
          return res;
        }

Solution graph:

.. image:: elasto-statics-sln.png

.. seealso::
  
   :ref:`example-heat-conduction`
