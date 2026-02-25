.. _changelog:


#########
Changelog
#########

******************
0.3.3 (unreleased)
******************

Internal Changes
================
- Removed the reduced_environment.txt file since it isn't needed anymore (:pull:`20`). By `Nathan Miller`_.

******************
0.3.2 (02-17-2025)
******************

Internal Changes
================
- Updated to use the conda-forge version of tardigrade_hydra (:pull:`18`). By `Nathan Miller`_.
- Updated to reflect upstream changes in tardigrade_hydra API (:pull:`19`). By `Nathan Miller`_.

Release
=======
- Release version (:pull:`20`). By `Nathan Miller`_.

******************
0.3.1 (02-13-2025)
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
