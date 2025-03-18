.. _changelog:


#########
Changelog
#########

******************
0.2.0 (11-07-2024)
******************

New Features
============
- Added iterator based calculation of the balance of mass (:pull:`9`). By `Nathan Miller`_.
- Added a multiphase balance of mass calculation (:pull:`10`). By `Nathan Miller`_.
- Added the calculation of the non-divergence part of the balance of linear momentum (:pull:`11`). By `Nathan Miller`_.
- Added the calculation of the divergence part of the balance of linear momentum (:pull:`11`). By `Nathan Miller`_.
- Added the multiphase balance of linear momentum (:pull:`12`). By `Nathan Miller`_.
- Added the calculation of the balance of energy for single and multiphase reacting continuum (:pull:`13`). By `Nathan Miller`_.
- Added a hexahedral finite element to test the balance equations in a FEA context (:pull:`15`). By `Nathan Miller`_.
- Added the FEA Jacobians for the balance of mass, linear momentum, and energy (:pull:`16`). By `Nathan Miller`_.
- Added the balance of volume fraction (:pull:`17`). By `Nathan Miller`_.
- Added material model version for the balance equations. (:pull:`18`). By `Nathan Miller`_.
- Added a constraint equation enabling complex relationships between the material state and the internal energy. (:pull:`19`). By `Nathan Miller`_.
- Added a constraint equation for the displacement. (:pull:`20`). by `Nathan Miller`_.
- Marked all equations as inline and corrected documentation build errors. (:pull:`21`). By `Nathan Miller`_.
- Changed the definition of the template parameter material_response_num_dof from being the total number of dof (which would change with phase count)
  to just being the number of dof per phase (including the additional dof). (:pull:`22`). By `Nathan Miller`_.

Internal Changes
================
- Corrected the nomenclature for the derivative of gradient quantities (:pull:`14`). By `Nathan Miller`_.
- Moved the declarations of the balance of linear momentum to the header file. (:pull:`23`). By `Nathan Miller`_.
- Corrected uninitialized warning in balance of mass. (:pull:`24`). By `Nathan Miller`_.
- Moved balance of volume fraction declaration to the header file (:pull:`26`). By `Nathan Miller`_.

Breaking Changes
================
- Changed the order of the templates for the balance of mass to be consistent with the other balance equations. (:pull:`25`). By `Nathan Miller`_.

******************
0.1.0 (11-07-2024)
******************

Release
=======
- Released version 0.1.0 (:pull:`8`). By `Nathan Miller`_.

New Features
============
- Initial commit of the balance equation repository (:pull:`1`). By `Nathan Miller`_.
- Removed mentions of tardigrade hydra from the readme (:pull:`2`). By `Nathan Miller`_.
- Added the calculation of the balance of mass (:pull:`3`). By `Nathan Miller`_.
- Added the calculation of the derivative of the spatial gradient of a quantity (:pull:`4`). By `Nathan Miller`_.
- Allow the version to be specified when doing a FetchContent build (:pull:`7`). By `Nathan Miller`_.

Internal Changes
================
- Removed shared only library output (:pull:`5`). By `Nathan Miller`_.

Bug Fixes
=========
- Removed leading whitespace for add_library call (:pull:`6`). By `Nathan Miller`_.
