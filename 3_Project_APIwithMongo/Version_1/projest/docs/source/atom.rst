atom
====

.. toctree::
   :maxdepth: 4

.. image:: ../images/atom.jpeg
    :width: 500px
    :align: center
    :height: 500px
    :alt: alternate text

information
-----------

.. automodule:: atom
    :members:
    :undoc-members:

* **atom_name:** line[12:16] from pdb file. You obtain all the name of the atom. C-α -> alfa carbon. String.
* **x_coordinate:** line[30:38] from pdb file. X coordinate. Float.
* **y_coordinate:** line[38:46] from pdb file. Y coordinate. Float.
* **z_coordinate:** line[46:54] from pdb file. Z coordinate. Float.
* **occupancy:** line[54:60] from pdb file. 1 corresponds to 100% occupancy. Float.
* **temperature_factor:** line[60:66] from pdb file. The temperature for each atom in the structure. Float.
* **element_symbol:** line[76:78] from pdb file. The element of the atom. String.
