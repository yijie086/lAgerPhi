// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
//
// This file is part of lAger.
//
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef LAGER_PHYSICS_FERMI_LOADED
#define LAGER_PHYSICS_FERMI_LOADED

// =============================================================================
// Fermi momentum distributions
// =============================================================================

namespace lager {
namespace physics {

// =============================================================================
// 1987 NBS Fermi-momentum implementation
// By J.S O'Connell and J.W. Lightbody, Jr
// See fermi87f for full attribution
double fermi87(double P, int A);

} // physics
} // lager

#endif
