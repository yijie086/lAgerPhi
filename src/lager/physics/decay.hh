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

#ifndef LAGER_PHYSICS_DECAY_LOADED
#define LAGER_PHYSICS_DECAY_LOADED

#include <lager/core/particle.hh>
#include <tuple>

namespace lager {
namespace physics {

// =============================================================================
// GENERIC TWO BODY DECAY
// =============================================================================

// two body decay of a particle 'part' into two particles xx (xx.first,
// xx.second), with angles of the first decay particle ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
void decay_2body(const particle& part, const double theta_1, const double phi_1,
                 std::pair<particle, particle>& xx);
// same, but will also store the decay particles in the helicity frame
void decay_2body(const particle& part, const double theta_1, const double phi_1,
                 std::pair<particle, particle>& xx,
                 std::pair<particle, particle>& xx_cm);

} // physics
} // lager

#endif
