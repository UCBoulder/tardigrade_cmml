.. _changelog:


#########
Changelog
#########

******************
0.3.7 (unreleased)
******************

******************
0.3.6 (04-01-2026)
******************

Bug Fixes
=========
- Fixed a bug in the linear test material where the response was not being initialized to zero (:pull:`28`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`29`). By `Nathan Miller`_.

******************
0.3.5 (03-27-2026)
******************

New Features
============
- Added a linear test material (:pull:`26`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`27`). By `Nathan Miller`_.

******************
0.3.4 (03-18-2026)
******************

New Features
============
- Added the CHIPFoam hyperelastic material model (:pull:`24`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`25`). By `Nathan Miller`_.

******************
0.3.3 (03-17-2026)
******************

Internal Changes
================
- Removed the reduced_environment.txt file since it isn't needed anymore (:pull:`20`). By `Nathan Miller`_.
- Removed reduced environment file (:pull:`21`). By `Nathan Miller`_.
- Updated to allow builds with Eigen3 >= 5 (:pull:`22`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`23`). By `Nathan Miller`_.

******************
0.3.2 (02-17-2026)
******************

Internal Changes
================
- Updated to use the conda-forge version of tardigrade_hydra (:pull:`18`). By `Nathan Miller`_.
- Updated to reflect upstream changes in tardigrade_hydra API (:pull:`19`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`20`). By `Nathan Miller`_.

******************
0.3.1 (02-13-2026)
******************

Internal Changes
================
- Incorporated breaking changes from tardigrade_hydra api (:pull:`15`). By `Nathan Miller`_.
- Updated for conda-packaging and using linters (:pull:`16`). By `Nathan Miller`_.
- Require tardigrade_hydra for the environment during testing (:pull:`18`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`17`). By `Nathan Miller`_.

******************
0.3.0 (09-17-2025)
******************

Bug Fixes
=========
- Fixed a bug in the the dof velocity driven deformation material model where the SDVs weren't being updated correctly (:pull:`9`). By `Nathan Miller`_.

New Features
============
- Added a mass diffusion model to DefinedPlasticEvolution (:pull:`11`). By `Nathan Miller`_.
- Enabled mass change rate and internal heat generation rate calculations (:pull:`12`). By `Nathan Miller`_.

Internal Changes
================
- Changed the default integrator for the defined deformation material model to fully implicit trapezoidal (:pull:`13`). By `Nathan Miller`_.
- Added additional outputs for the testing (:pull:`13`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`14`). By `Nathan Miller`_.

******************
0.2.0 (06-12-2025)
******************

New Features
============
- Added a dof velocity gradient driven deformation material model DefinedPlasticEvolution (:pull:`7`). By `Nathan Miller`_.

Breaking Changes
================
- Changed BasicReactingSolid to BasicSolid since it is more accurate (:pull:`7`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`8`). By `Nathan Miller`_.

******************
0.1.0 (05-21-2025)
******************

New Features
============
- Initial commit (:pull:`1`). By `Nathan Miller`_.
- Added the basic reacting continuum class (:pull:`2`). By `Nathan Miller`_.
- Added the full material model library build (:pull:`3`). By `Nathan Miller`_.
- Added a function that returns the evaluate model results size (:pull:`4`). By `Nathan Miller`_.

Internal Changes
================
- Prepared the library for conda packaging (:pull:`5`). By `Nathan Miller`_.


Release
=======
- Release version (:pull:`6`). By `Nathan Miller`_.
