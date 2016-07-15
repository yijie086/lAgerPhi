#ifndef PCSIM_PHYSICS_CONSTANTS_LOADED
#define PCSIM_PHYSICS_CONSTANTS_LOADED

namespace pcsim {
namespace physics {

// Fine structure
constexpr const double ALPHA = 0.0072973525664;

// Mass (in [GeV])
constexpr const double M_ELECTRON = 0.000510999; // electron mass
constexpr const double M_PROTON = 0.938272;      // proton mass
constexpr const double M_JPSI = 3.09692;         // J/Psi pole mass

// Mass squared
constexpr const double M2_ELECTRON = M_ELECTRON * M_ELECTRON; // electron
constexpr const double M2_PROTON = M_PROTON * M_PROTON;       // proton
constexpr const double M2_JPSI = M_JPSI * M_JPSI;             // J/Psi

// 1/(2Mass)
constexpr const double ONE_OVER_2M_PROTON = 1. / (2. * M_PROTON); // proton

// Width
constexpr const double W_JPSI = 92.9e-6; // J/Psi width

// Branching Ratio
constexpr const double B_JPSI_ELEC = 0.0602; // J/Psi->e+e- branching ratio

} // pcsim
} // physics

#endif
