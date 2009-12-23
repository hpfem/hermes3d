.. _example-heat-conduction:

Heat Conduction
===============

.. index::
   single: mesh; fixed
   single: problem; time-depenedent
   single: problem; elliptic

This example shows how to implement a solver for time dependent problem on a fixed mesh.

PDE solved:

.. math::
   :nowrap:

   \begin{eqnarray*}
   \frac{\partial u}{\partial t} - \Delta u &= f &\hbox{ in }\Omega \\
                                          u &= 0 &\hbox{ on }\partial\Omega
   \end{eqnarray*}

Exact solution is:

.. math:: u = \sin(t) (1 - x^2) (1 - y^2) (1 - z^2)


.. literalinclude:: ../../examples/heat-conduction/main.cc
   :language: c
   :linenos:
   :lines: 22-

.. seealso::
  
   :ref:`example-fichera-corner`
