.. _example-bessel:

Bessel
======

.. index::
   single: mesh; fixed
   single: problem; elliptic
   single: Hcurl


This example shows how to solve complex valued problem.

PDE solved:

.. math::
   :nowrap:

   \begin{eqnarray*}
   \nabla \times (\nabla \times {\mathbf E}) - {\mathbf E} &= {\mathbf F} &\hbox{ in }\Omega \\
            {\mathbf E} \times \nu &= 0 &\hbox{ on }\Gamma_P\\
            (\nabla \times {\mathbf E}) \times \nu - i {\mathbf E} &= g &\hbox{ on }\Gamma_I
   \end{eqnarray*}

Exact solution is:

.. math:: {\mathbf E} = \nabla \times (J_{\frac{2}{3}}(r) \cos(\frac{2}{3} \theta))


Solution:

.. image:: bessel-sln.png


Source code:

.. literalinclude:: ../../examples/bessel/main.cc
   :language: c
   :linenos:
   :lines: 23-
