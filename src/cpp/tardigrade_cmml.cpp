/*!
=====================================================================
|                        tardigrade_cmml.cpp                        |
=====================================================================
| (C)ollected (M)aterial (M)odel (L)ibrary                          |
| prounced Camel                                                    |
|                                                                   |
| A header file which defines a class which registers all of the    |
| available material models and their interface. The library is     |
| implemented in such a way that someone can register a new         |
| constitutive model and have it available for use merely by        |
| calling it by name. This allows us to re-use code as much as      |
| possible without needing to do time consuming rewrites of already |
| existing code.                                                    |
---------------------------------------------------------------------
| Note: Registration approach taken from stackoverflow question     |
|       compile time plugin system 2                                |
=====================================================================
*/

#include "tardigrade_cmml.h"
