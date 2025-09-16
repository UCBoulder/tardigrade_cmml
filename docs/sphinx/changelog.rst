.. _changelog:


#########
Changelog
#########

******************
0.2.1 (unreleased)
******************

Bug Fixes
=========
- Fixed a bug in the the dof velocity driven deformation material model where the SDVs weren't being updated correctly (:pull:`9`). By `Nathan Miller`_.

New Features
============
- Added a mass diffusion model to DefinedPlasticEvolution (:pull:`11`). By `Nathan Miller`_.
- Enabled mass change rate and internal heat generation rate calculations (:pull:`12`). By `Nathan Miller`_.

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
