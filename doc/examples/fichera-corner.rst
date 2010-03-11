.. _example-fichera-corner:

Fichera Corner
==============

.. index::
   single: mesh; dynamical
   single: problem; elliptic


This example shows how to implement the adaptivity loop.

PDE solved:

.. math::
   :nowrap:

   \begin{eqnarray*}
   - \Delta u &= f &\hbox{ in }\Omega \\
            u &= g &\hbox{ on }\partial\Omega
   \end{eqnarray*}

Exact solution is:

.. math:: u = (x^2 + y^2 + z^2)^{0.25}


Convergence graphs:

.. image:: fichera-conv.png

.. image:: fichera-conv-time.png


Solution and hp-mesh:

.. image:: fichera-sln.png

.. image:: fichera-order.png


Source code:

.. literalinclude:: ../../examples/fichera/main.cc
   :language: c
   :linenos:
   :lines: 30-

.. seealso::

   :ref:`example-heat-conduction`
