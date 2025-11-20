#!/usr/bin/env python3
import uproot
import numpy as np

def main():
    # --- CONFIGURATION ---
    fname  = "outputs/hatta-ep-phi-10GeV.ep-phi.4pi.run00484-lumi100.root"
    tree_name = "lAger"
    Ebeam  = 10.6            # GeV
    Mp     = 0.938272        # proton mass [GeV]
    alpha  = 1.0/137.0
    L      = 1e17            # integrated luminosity [b^-1] (1000 fb^-1)

    # bin edges
    Q2_lo, Q2_hi = 3.0, 4.0
    W_lo,  W_hi  = 2.2, 2.3
    t_lo,  t_hi  = -1.1, -0.8

    # bin widths
    dQ2 = Q2_hi - Q2_lo
    dW  = W_hi  - W_lo
    dt  = t_hi  - t_lo

    # --- READ TTree ---
    f = uproot.open(fname)
    tree = f[tree_name]
    # pull only the branches we need
    arr = tree.arrays(["Q2","W","t","nu","epsilon"], library="np")

    # --- SELECT BIN ---
    mask = (
        (arr["Q2"] > Q2_lo) & (arr["Q2"] < Q2_hi) &
        (arr["W"]  >  W_lo) & (arr["W"]  <  W_hi ) &
        (arr["t"]  >  t_lo ) & (arr["t"]  <  t_hi )
    )
    N = mask.sum()
    if N == 0:
        print("No events in that bin — check your boundaries.")
        return

    Q2      = arr["Q2"][mask]
    W       = arr["W" ][mask]
    nu      = arr["nu"][mask]
    eps     = arr["epsilon"][mask]

    # --- VIRTUAL-PHOTON FLUX (Hand convention) ---
    # Kinematic factors event-by-event
    Eprime = Ebeam - nu
    K      = (W**2 - Mp**2) / (2*Mp)
    # Γ_v = α/(2π^2) * (K / Q2) * (E'/E) * 1/(1-ε)
    Gamma  = alpha / (2*np.pi**2) * (K/Q2) * (Eprime/Ebeam) * (1.0/(1.0 - eps))
    Gavg   = np.mean(Gamma)

    # --- REDUCED CROSS SECTION σ_red ---
    # σ_red = N / (L * ΔW * ΔQ2 * Δt * <Γ_v>)
    sigma_red_barn = N / (L * dW * dQ2 * dt * Gavg)
    sigma_red_nb   = sigma_red_barn * 1e9   # 1 b = 1e9 nb

    # --- OUTPUT ---
    print(f"Events in bin         : N = {N}")
    print(f"Bin widths            : ΔW={dW:.3f}, ΔQ2={dQ2:.3f}, Δt={dt:.3f}")
    print(f"Average flux <Γ_v>    : {Gavg:.3e} GeV⁻²")
    print(f"Reduced σ (barn)      : {sigma_red_barn:.3e} b")
    print(f"Reduced σ (nanobarn)  : {sigma_red_nb:.3e} nb")

if __name__ == "__main__":
    main()
