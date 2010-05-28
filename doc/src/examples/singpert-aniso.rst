.. _example-sing-pert:

Singular Perturbation
=====================

.. index::
   single: mesh; dynamical
   single: problem; elliptic


This examples demonstrates how the anisotropic refinements can save a big amount of degrees
of freedom.


PDE solved:

.. math::
   :nowrap:

   \begin{eqnarray*}
   - \Delta u + K^2 u &= F &\hbox{ in }\Omega \\
                    u &= 0 &\hbox{ on }\partial\Omega
   \end{eqnarray*}

Convergence graphs:

.. image:: singpert-aniso-conv.png

.. image:: singpert-aniso-conv-time.png


Solution and hp-mesh:

.. image:: singpert-aniso-sln.png

.. image:: singpert-aniso-order.png




.. literalinclude:: ../../examples/singpert-aniso/main.cc
   :language: c
   :linenos:
   :lines: 30-

.. seealso::
  
   :ref:`example-fichera-corner`
