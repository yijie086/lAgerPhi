// lundANA.C
// Read LUND file with ep -> ep K+ K- events,
// store beam energy and e-, p, K+, K- four-momenta and derived kinematics
// (phi_deg, theta_deg, p, missing masses, Emiss, PTmiss, Q2, xB, t, W, m_KK, phi_Trento)
// into a TTree.
//
// Angles are stored in DEGREES.

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Compute Trento phi angle (in degrees) between leptonic (k, k') and hadronic (q, V) planes.
// Definition: see Bacchetta et al., Trento conventions.
// Input 4-vectors are in lab; only spatial parts are used.
double computePhiTrentoDeg(const TLorentzVector &lv_k,
                           const TLorentzVector &lv_kp,
                           const TLorentzVector &lv_q,
                           const TLorentzVector &lv_V)
{
    TVector3 k  = lv_k.Vect();
    TVector3 kp = lv_kp.Vect();
    TVector3 q  = lv_q.Vect();
    TVector3 V  = lv_V.Vect();

    // z-axis opposite to q direction
    TVector3 z = -q.Unit();

    TVector3 n_l = k.Cross(kp);   // normal to leptonic plane
    TVector3 n_h = q.Cross(V);    // normal to hadronic plane

    if (n_l.Mag() == 0 || n_h.Mag() == 0) {
        return -999.0; // undefined
    }

    n_l = n_l.Unit();
    n_h = n_h.Unit();

    double cosphi = n_l.Dot(n_h);
    if (cosphi >  1.0) cosphi =  1.0;
    if (cosphi < -1.0) cosphi = -1.0;

    double phi = std::acos(cosphi); // in [0, pi]

    // Sign using (n_l x n_h) · z
    double sign = z.Dot(n_l.Cross(n_h));
    if (sign < 0) {
        phi = 2.0 * M_PI - phi;
    }

    return phi * 180.0 / M_PI; // degrees in [0,360)
}

void lundANA(const char *inFile  = "CLAS-ep-phi-10GeV.ep-phi.4pi.run00024-5000_clean.txt",
             const char *outFile = "events.root")
{
    std::ifstream fin(inFile);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open input LUND file: " << inFile << std::endl;
        return;
    }

    TFile *fout = new TFile(outFile, "RECREATE");
    if (fout->IsZombie()) {
        std::cerr << "Error: cannot create output ROOT file: " << outFile << std::endl;
        return;
    }

    // Beam energy (from header)
    double beamE = 0.0;

    // 4-momenta for e-, p, K+, K-
    double e_px, e_py, e_pz, e_E;
    double p_px, p_py, p_pz, p_E;
    double kp_px, kp_py, kp_pz, kp_E;
    double km_px, km_py, km_pz, km_E;

    // Derived single-particle kinematics (angle in DEG)
    double e_p,  e_theta_deg,  e_phi_deg;
    double p_p,  p_theta_deg,  p_phi_deg;
    double kp_p, kp_theta_deg, kp_phi_deg;
    double km_p, km_theta_deg, km_phi_deg;

    // Missing mass^2
    double mm2_ep;
    double mm2_epKpKm;
    double mm2_eKpKm;
    double mm2_epKp;
    double mm2_epKm;

    // Missing energy and transverse missing momentum
    double E_miss;
    double PT_miss;

    // DIS variables and KK invariant mass
    double Q2;
    double xB;
    double t;
    double W;
    double m_KK;

    // Trento phi for K+K- system (degrees)
    double phi_trento_KK_deg;

    // Constants: electron and proton mass in GeV
    const double me = 0.00051099895;
    const double mp = 0.938272088;

    // Create tree and branches
    TTree *T = new TTree("T", "ep -> ep K+ K- generated kinematics");

    // Beam
    T->Branch("beamE", &beamE, "beamE/D");

    // Raw 4-momenta
    T->Branch("e_px",  &e_px,  "e_px/D");
    T->Branch("e_py",  &e_py,  "e_py/D");
    T->Branch("e_pz",  &e_pz,  "e_pz/D");
    T->Branch("e_E",   &e_E,   "e_E/D");

    T->Branch("p_px",  &p_px,  "p_px/D");
    T->Branch("p_py",  &p_py,  "p_py/D");
    T->Branch("p_pz",  &p_pz,  "p_pz/D");
    T->Branch("p_E",   &p_E,   "p_E/D");

    T->Branch("kp_px", &kp_px, "kp_px/D");
    T->Branch("kp_py", &kp_py, "kp_py/D");
    T->Branch("kp_pz", &kp_pz, "kp_pz/D");
    T->Branch("kp_E",  &kp_E,  "kp_E/D");

    T->Branch("km_px", &km_px, "km_px/D");
    T->Branch("km_py", &km_py, "km_py/D");
    T->Branch("km_pz", &km_pz, "km_pz/D");
    T->Branch("km_E",  &km_E,  "km_E/D");

    // Single-particle kinematics (DEGREES)
    T->Branch("e_p",         &e_p,         "e_p/D");
    T->Branch("e_theta_deg", &e_theta_deg, "e_theta_deg/D");
    T->Branch("e_phi_deg",   &e_phi_deg,   "e_phi_deg/D");

    T->Branch("p_p",         &p_p,         "p_p/D");
    T->Branch("p_theta_deg", &p_theta_deg, "p_theta_deg/D");
    T->Branch("p_phi_deg",   &p_phi_deg,   "p_phi_deg/D");

    T->Branch("kp_p",         &kp_p,         "kp_p/D");
    T->Branch("kp_theta_deg", &kp_theta_deg, "kp_theta_deg/D");
    T->Branch("kp_phi_deg",   &kp_phi_deg,   "kp_phi_deg/D");

    T->Branch("km_p",         &km_p,         "km_p/D");
    T->Branch("km_theta_deg", &km_theta_deg, "km_theta_deg/D");
    T->Branch("km_phi_deg",   &km_phi_deg,   "km_phi_deg/D");

    // Missing mass^2
    T->Branch("mm2_ep",     &mm2_ep,     "mm2_ep/D");
    T->Branch("mm2_epKpKm", &mm2_epKpKm, "mm2_epKpKm/D");
    T->Branch("mm2_eKpKm",  &mm2_eKpKm,  "mm2_eKpKm/D");
    T->Branch("mm2_epKp",   &mm2_epKp,   "mm2_epKp/D");
    T->Branch("mm2_epKm",   &mm2_epKm,   "mm2_epKm/D");

    // Missing energy and PT
    T->Branch("E_miss",  &E_miss,  "E_miss/D");
    T->Branch("PT_miss", &PT_miss, "PT_miss/D");

    // DIS variables and KK invariant mass
    T->Branch("Q2",   &Q2,   "Q2/D");
    T->Branch("xB",   &xB,   "xB/D");
    T->Branch("t",    &t,    "t/D");
    T->Branch("W",    &W,    "W/D");
    T->Branch("m_KK", &m_KK, "m_KK/D");

    // Trento phi (degrees)
    T->Branch("phi_trento_KK_deg", &phi_trento_KK_deg, "phi_trento_KK_deg/D");

    std::string line;
    Long64_t nEvents = 0;
    Long64_t nFilled = 0;

    // Main loop over events
    while (true) {
        if (!std::getline(fin, line)) break;   // EOF
        if (line.empty() || line[0] == '#') continue;

        std::istringstream hs(line);

        int    npart;
        int    h_i2, h_i3, h_beamPid, h_i8, h_i9;
        double h_d4, h_d5, h_weight;

        hs >> npart >> h_i2 >> h_i3
           >> h_d4 >> h_d5
           >> h_beamPid >> beamE
           >> h_i8 >> h_i9 >> h_weight;

        if (!hs) continue;
        ++nEvents;

        // Zero out everything
        e_px = e_py = e_pz = e_E = 0.0;
        p_px = p_py = p_pz = p_E = 0.0;
        kp_px = kp_py = kp_pz = kp_E = 0.0;
        km_px = km_py = km_pz = km_E = 0.0;

        bool has_e=false, has_p=false, has_kp=false, has_km=false;

        // Read particle lines
        for (int i = 0; i < npart; ++i) {
            if (!std::getline(fin, line)) break;
            if (line.empty() || line[0]=='#') { --i; continue; }

            std::istringstream ps(line);

            int idx, ist, pid, parent, daughter;
            double dummy2;
            double px, py, pz, E, m, vx, vy, vz;

            ps >> idx >> dummy2 >> ist >> pid >> parent >> daughter
               >> px >> py >> pz >> E >> m >> vx >> vy >> vz;

            if (!ps) continue;

            if (pid == 11) {
                e_px = px; e_py = py; e_pz = pz; e_E = E; has_e = true;
            } else if (pid == 2212) {
                p_px = px; p_py = py; p_pz = pz; p_E = E; has_p = true;
            } else if (pid == 321) {
                kp_px = px; kp_py = py; kp_pz = pz; kp_E = E; has_kp = true;
            } else if (pid == -321) {
                km_px = px; km_py = py; km_pz = pz; km_E = E; has_km = true;
            }
        }

        if (!(has_e && has_p && has_kp && has_km)) continue;

        // Build TLorentzVectors
        double p_beam_mag = std::sqrt(beamE*beamE - me*me);
        TLorentzVector lv_beam(0, 0, p_beam_mag, beamE);
        TLorentzVector lv_target(0, 0, 0, mp);

        TLorentzVector lv_e(e_px, e_py, e_pz, e_E);
        TLorentzVector lv_p(p_px, p_py, p_pz, p_E);
        TLorentzVector lv_kp(kp_px, kp_py, kp_pz, kp_E);
        TLorentzVector lv_km(km_px, km_py, km_pz, km_E);

        // Single-particle kinematics (DEGREES)
        auto fill_deg = [](const TLorentzVector &lv,
                           double &p,
                           double &theta_deg,
                           double &phi_deg)
        {
            p         = lv.P();
            theta_deg = lv.Theta() * 180.0 / M_PI;
            phi_deg   = lv.Phi()   * 180.0 / M_PI;
        };

        fill_deg(lv_e,  e_p,  e_theta_deg,  e_phi_deg);
        fill_deg(lv_p,  p_p,  p_theta_deg,  p_phi_deg);
        fill_deg(lv_kp, kp_p, kp_theta_deg, kp_phi_deg);
        fill_deg(lv_km, km_p, km_theta_deg, km_phi_deg);

        // Missing masses
        TLorentzVector lv_init = lv_beam + lv_target;
        mm2_ep     = (lv_init - lv_e - lv_p).M2();
        mm2_epKpKm = (lv_init - lv_e - lv_p - lv_kp - lv_km).M2();
        mm2_eKpKm  = (lv_init - lv_e - lv_kp - lv_km).M2();
        mm2_epKp   = (lv_init - lv_e - lv_p - lv_kp).M2();
        mm2_epKm   = (lv_init - lv_e - lv_p - lv_km).M2();

        // Missing energy and PT
        TLorentzVector lv_miss = lv_init - lv_e - lv_p - lv_kp - lv_km;
        E_miss  = lv_miss.E();
        PT_miss = std::sqrt(lv_miss.Px()*lv_miss.Px() +
                            lv_miss.Py()*lv_miss.Py());

        // DIS variables
        TLorentzVector lv_q = lv_beam - lv_e;        // virtual photon
        Q2 = -lv_q.M2();

        double nu = beamE - e_E;
        if (nu > 0.0) {
            xB = Q2 / (2.0 * mp * nu);
        } else {
            xB = -999.0;
        }

        TLorentzVector lv_W = lv_target + lv_q;
        W = lv_W.M();

        TLorentzVector lv_t = lv_p - lv_target;
        t = lv_t.M2();

        TLorentzVector lv_KK = lv_kp + lv_km;
        m_KK = lv_KK.M();

        // Trento phi for KK system
        phi_trento_KK_deg = computePhiTrentoDeg(lv_beam, lv_e, lv_q, lv_KK);

        T->Fill();
        ++nFilled;
    }

    fout->cd();
    T->Write();
    fout->Close();
    fin.close();

    std::cout << "Finished. Headers read: " << nEvents
              << ", entries filled: " << nFilled << std::endl;
}
