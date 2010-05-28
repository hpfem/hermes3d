.. _example-layer:

Layer
=====

.. index::
   single: mesh; dynamical
   single: problem; elliptic


PDE solved:

.. math::
   :nowrap:

   \begin{eqnarray*}
   - \Delta u &= f &\hbox{ in }\Omega \\
            u &= g &\hbox{ on }\partial\Omega
   \end{eqnarray*}

Exact solution is:

.. math:: u = atan(k \cdot \sqrt{(x + 0.25)^2 + (y + 0.25)^2 + (z + 0.25)^2} - \frac{\pi}{3})

Convergence graphs:

.. image:: layer-conv.png

.. image:: layer-conv-time.png


Solution and hp-mesh:

.. image:: layer-sln.png

.. image:: layer-order.png


.. literalinclude:: ../../examples/layer/main.cc
   :language: c
   :linenos:
   :lines: 30-

.. seealso::

   :ref:`example-fichera-corner`

