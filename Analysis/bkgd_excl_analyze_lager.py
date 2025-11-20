import ROOT
import numpy as np
import sys

def main(input_file, output_file):
    # Open the ROOT file
    file = ROOT.TFile(input_file)

    # Get the TTree
    tree = file.Get("DiffPhiBkgdEvents")

    # Initialize lists to store momenta and polar angles
    electron_momenta = []
    proton_momenta = []
    kaon_plus_momenta = []
    kaon_minus_momenta = []
    
    electron_theta = []
    proton_theta = []
    kaon_plus_theta = []
    kaon_minus_theta = []

    # Kinematic variables
    W_values = []
    Q2_values = []
    x_values = []
    t_values = []

    filtered_electron_momenta = []
    filtered_proton_momenta = []
    filtered_kaon_plus_momenta = []
    filtered_kaon_minus_momenta = []

    filtered_electron_theta = []
    filtered_proton_theta = []
    filtered_kaon_plus_theta = []
    filtered_kaon_minus_theta = []

    # Kinematic variables after filtering
    filtered_W_values = []
    filtered_Q2_values = []
    filtered_x_values = []
    filtered_t_values = []

    # Resolutions
    Q2_resolutions = []
    W_resolutions = []
    x_resolutions = []
    y_resolutions = []
    t_resolutions = []

    # Missing mass
    missing_mass_values = []

    KK_mass_values = []

    # Cuts (modify these values as needed)
    W_cut = (1.0, 2.4)
    Q2_cut = (3.5, 60.0)
    x_cut = (0.0001, 0.9999)

    # Define initial state four-momentum (proton at rest, incoming electron)
    initial_proton = ROOT.TLorentzVector(0, 0, 0, 0.938)  # Proton mass in GeV
    initial_electron = ROOT.TLorentzVector(0, 0, 10.6, np.sqrt(10.6**2 + 0.000511**2))  # Incoming electron 10.6 GeV in z-direction

    # Function to calculate polar angle theta in degrees
    def calculate_theta(p):
        return np.degrees(np.arccos(p.Pz() / p.P()))

    # Function to calculate kinematic variables
    def calculate_kinematics(electron):
        q = initial_electron - electron
        Q2 = -q.M2()
        W = (initial_proton + q).M()
        x = Q2 / (2 * initial_proton.Dot(q))
        y = Q2 / (2 * initial_electron.Dot(q))
        return W, Q2, x, y

    # Smearing functions
    def smear_momentum(momentum, resolution=0.02):
        return np.random.normal(momentum, resolution * momentum)

    def smear_angle(angle, resolution=0.001):
        return angle + np.random.normal(0, resolution)

    # Open the detector acceptance ROOT file
    acceptance_file_path = "/home/klest/pcsim_proc_detector_solid-fastmc_acceptance_solid_JPsi_electron_target315_output.root"
    acceptance_file = ROOT.TFile(acceptance_file_path)
    overall_acceptance_hist = acceptance_file.Get("acceptance_ThetaP_overall")
    charged_acceptance_hist = acceptance_file.Get("acceptance_ThetaP_forwardangle")

    def scale_histogram(hist, max_value):
        hist_max = hist.GetMaximum()
        scale_factor = max_value / hist_max
        hist.Scale(scale_factor)

    # Scale the acceptance histograms
    scale_histogram(overall_acceptance_hist, 0.9)
    scale_histogram(charged_acceptance_hist, 0.9)
    iteration_count = 0

    # Loop over the entries in the tree
    for event in tree:
        electron = event.Electron
        charged_particles = [
            (event.ChargedParticle1, event.ChargedPID1),
            (event.ChargedParticle2, event.ChargedPID2),
            (event.ChargedParticle3, event.ChargedPID3),
            (event.ChargedParticle4, event.ChargedPID4),
            (event.ChargedParticle5, event.ChargedPID5),
        ]

        proton = None
        kaon_plus = None
        kaon_minus = None

        for charged_particle, pid in charged_particles:
            if pid == 2212:  # Proton
                proton = charged_particle
            elif pid == -211:  # K+
                kaon_plus = charged_particle
            elif pid == 321:  # K-
                kaon_minus = charged_particle

        if not proton or not kaon_plus or not kaon_minus:
            continue

        # Calculate kinematic variables
        W, Q2, x, y = calculate_kinematics(electron)
        t = (initial_proton - proton).M2()

        # Apply generator-level cuts
        if (W_cut[0] <= W <= W_cut[1]) and (Q2_cut[0] <= Q2 <= Q2_cut[1]) and (x_cut[0] <= x <= x_cut[1]):
            # Append generator-level kinematic variables
            W_values.append(W)
            Q2_values.append(Q2)
            x_values.append(x)
            t_values.append(t)

            # Append generator-level momenta and polar angles for all particles
            electron_momenta.append(electron.P())
            proton_momenta.append(proton.P())
            kaon_plus_momenta.append(kaon_plus.P())
            kaon_minus_momenta.append(kaon_minus.P())

            electron_theta.append(calculate_theta(electron))
            proton_theta.append(calculate_theta(proton))
            kaon_plus_theta.append(calculate_theta(kaon_plus))
            kaon_minus_theta.append(calculate_theta(kaon_minus))

        # Apply acceptance cuts
        acceptances = []

        # Electron acceptance
        electron_theta_val = calculate_theta(electron)
        electron_momentum_val = electron.P()
        electron_bin_x = overall_acceptance_hist.GetXaxis().FindBin(electron_theta_val)
        electron_bin_y = overall_acceptance_hist.GetYaxis().FindBin(electron_momentum_val)
        electron_acceptance = overall_acceptance_hist.GetBinContent(electron_bin_x, electron_bin_y)
        if np.random.random() < electron_acceptance:
            acceptances.append(True)
        else:
            acceptances.append(False)

        # Proton acceptance
        proton_theta_val = calculate_theta(proton)
        proton_momentum_val = proton.P()
        proton_bin_x = charged_acceptance_hist.GetXaxis().FindBin(proton_theta_val)
        proton_bin_y = charged_acceptance_hist.GetYaxis().FindBin(proton_momentum_val)
        proton_acceptance = charged_acceptance_hist.GetBinContent(proton_bin_x, proton_bin_y)
        if np.random.random() < proton_acceptance:
            acceptances.append(True)
        else:
            acceptances.append(False)

        # Kaon+ acceptance
        kaon_plus_theta_val = calculate_theta(kaon_plus)
        kaon_plus_momentum_val = kaon_plus.P()
        kaon_plus_bin_x = charged_acceptance_hist.GetXaxis().FindBin(kaon_plus_theta_val)
        kaon_plus_bin_y = charged_acceptance_hist.GetYaxis().FindBin(kaon_plus_momentum_val)
        kaon_plus_acceptance = charged_acceptance_hist.GetBinContent(kaon_plus_bin_x, kaon_plus_bin_y)
        if np.random.random() < kaon_plus_acceptance:
            acceptances.append(True)
        else:
            acceptances.append(False)

        # Kaon- acceptance
        kaon_minus_theta_val = calculate_theta(kaon_minus)
        kaon_minus_momentum_val = kaon_minus.P()
        kaon_minus_bin_x = charged_acceptance_hist.GetXaxis().FindBin(kaon_minus_theta_val)
        kaon_minus_bin_y = charged_acceptance_hist.GetYaxis().FindBin(kaon_minus_momentum_val)
        kaon_minus_acceptance = charged_acceptance_hist.GetBinContent(kaon_minus_bin_x, kaon_minus_bin_y)
        if np.random.random() < kaon_minus_acceptance:
            acceptances.append(True)
        else:
            acceptances.append(False)

        # Event passes if electron, proton, kaon+ and kaon- are reconstructed
        if all(acceptances):
            smeared_electron_momentum = smear_momentum(electron.P())
            smeared_electron_theta = smear_angle(np.radians(calculate_theta(electron)))
            smeared_electron_phi = smear_angle(electron.Phi())

            smeared_electron_px = smeared_electron_momentum * np.sin(smeared_electron_theta) * np.cos(smeared_electron_phi)
            smeared_electron_py = smeared_electron_momentum * np.sin(smeared_electron_theta) * np.sin(smeared_electron_phi)
            smeared_electron_pz = smeared_electron_momentum * np.cos(smeared_electron_theta)
            smeared_electron_E = np.sqrt(smeared_electron_px**2 + smeared_electron_py**2 + smeared_electron_pz**2 + 0.000511**2)
            smeared_electron = ROOT.TLorentzVector()
            smeared_electron.SetPxPyPzE(smeared_electron_px, smeared_electron_py, smeared_electron_pz, smeared_electron_E)

            smeared_proton_momentum = smear_momentum(proton.P())
            smeared_proton_theta = smear_angle(np.radians(calculate_theta(proton)))
            smeared_proton_phi = smear_angle(proton.Phi())

            smeared_proton_px = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.cos(smeared_proton_phi)
            smeared_proton_py = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.sin(smeared_proton_phi)
            smeared_proton_pz = smeared_proton_momentum * np.cos(smeared_proton_theta)
            smeared_proton_E = np.sqrt(smeared_proton_px**2 + smeared_proton_py**2 + smeared_proton_pz**2 + 0.938**2)
            smeared_proton = ROOT.TLorentzVector()
            smeared_proton.SetPxPyPzE(smeared_proton_px, smeared_proton_py, smeared_proton_pz, smeared_proton_E)

            smeared_kaon_plus_momentum = smear_momentum(kaon_plus.P())
            smeared_kaon_plus_theta = smear_angle(np.radians(calculate_theta(kaon_plus)))
            smeared_kaon_plus_phi = smear_angle(kaon_plus.Phi())

            smeared_kaon_plus_px = smeared_kaon_plus_momentum * np.sin(smeared_kaon_plus_theta) * np.cos(smeared_kaon_plus_phi)
            smeared_kaon_plus_py = smeared_kaon_plus_momentum * np.sin(smeared_kaon_plus_theta) * np.sin(smeared_kaon_plus_phi)
            smeared_kaon_plus_pz = smeared_kaon_plus_momentum * np.cos(smeared_kaon_plus_theta)
            smeared_kaon_plus_E = np.sqrt(smeared_kaon_plus_px**2 + smeared_kaon_plus_py**2 + smeared_kaon_plus_pz**2 + 0.4937**2)  # K+ mass
            smeared_kaon_plus = ROOT.TLorentzVector()
            smeared_kaon_plus.SetPxPyPzE(smeared_kaon_plus_px, smeared_kaon_plus_py, smeared_kaon_plus_pz, smeared_kaon_plus_E)

            smeared_kaon_minus_momentum = smear_momentum(kaon_minus.P())
            smeared_kaon_minus_theta = smear_angle(np.radians(calculate_theta(kaon_minus)))
            smeared_kaon_minus_phi = smear_angle(kaon_minus.Phi())

            smeared_kaon_minus_px = smeared_kaon_minus_momentum * np.sin(smeared_kaon_minus_theta) * np.cos(smeared_kaon_minus_phi)
            smeared_kaon_minus_py = smeared_kaon_minus_momentum * np.sin(smeared_kaon_minus_theta) * np.sin(smeared_kaon_minus_phi)
            smeared_kaon_minus_pz = smeared_kaon_minus_momentum * np.cos(smeared_kaon_minus_theta)
            smeared_kaon_minus_E = np.sqrt(smeared_kaon_minus_px**2 + smeared_kaon_minus_py**2 + smeared_kaon_minus_pz**2 + 0.4937**2)  # K- mass
            smeared_kaon_minus = ROOT.TLorentzVector()
            smeared_kaon_minus.SetPxPyPzE(smeared_kaon_minus_px, smeared_kaon_minus_py, smeared_kaon_minus_pz, smeared_kaon_minus_E)

            filtered_electron_momenta.append(smeared_electron.P())
            filtered_electron_theta.append(np.degrees(smeared_electron_theta))
            filtered_W, filtered_Q2, filtered_x, filtered_y = calculate_kinematics(smeared_electron)
            
            if (W_cut[0] <= filtered_W <= W_cut[1]) and (Q2_cut[0] <= filtered_Q2 <= Q2_cut[1]) and (x_cut[0] <= filtered_x <= x_cut[1]):
                filtered_W_values.append(filtered_W)
                filtered_Q2_values.append(filtered_Q2)
                filtered_x_values.append(filtered_x)
                
                # Calculate resolutions
                Q2_resolution = filtered_Q2 - Q2
                W_resolution = filtered_W - W
                x_resolution = filtered_x - x
                y_resolution = filtered_y - y
                t_resolution = (initial_proton - smeared_proton).M2() - t
                
                Q2_resolutions.append(Q2_resolution)
                W_resolutions.append(W_resolution)
                x_resolutions.append(x_resolution)
                y_resolutions.append(y_resolution)
                t_resolutions.append(t_resolution)
                
                filtered_proton_momenta.append(smeared_proton.P())
                filtered_proton_theta.append(np.degrees(smeared_proton_theta))
                
                filtered_kaon_plus_momenta.append(smeared_kaon_plus.P())
                filtered_kaon_plus_theta.append(np.degrees(smeared_kaon_plus_theta))
                
                filtered_kaon_minus_momenta.append(smeared_kaon_minus.P())
                filtered_kaon_minus_theta.append(np.degrees(smeared_kaon_minus_theta))
                
                KK_mass = (smeared_kaon_plus+smeared_kaon_minus).M()
                KK_mass_values.append(KK_mass)
                # Calculate missing mass
                reconstructed_particles = [smeared_electron, smeared_proton, smeared_kaon_plus, smeared_kaon_minus]
                missing_mass = ((initial_electron + initial_proton) - sum(reconstructed_particles, ROOT.TLorentzVector())).M()
                missing_mass_values.append(missing_mass)

    # Save histograms to ROOT file
    output_root_file = ROOT.TFile(output_file, "RECREATE")

    hist_electron_momenta = ROOT.TH1F("electron_momenta", "Electron Momenta", 500, 0, 10.6)
    hist_proton_momenta = ROOT.TH1F("proton_momenta", "Proton Momenta", 500, 0, 10.6)
    hist_kaon_plus_momenta = ROOT.TH1F("kaon_plus_momenta", "Kaon+ Momenta", 500, 0, 10.6)
    hist_kaon_minus_momenta = ROOT.TH1F("kaon_minus_momenta", "Kaon- Momenta", 500, 0, 10.6)
    hist_electron_theta = ROOT.TH1F("electron_theta", "Electron Theta", 500, 0, 50)
    hist_proton_theta = ROOT.TH1F("proton_theta", "Proton Theta", 500, 0, 50)
    hist_kaon_plus_theta = ROOT.TH1F("kaon_plus_theta", "Kaon+ Theta", 500, 0, 50)
    hist_kaon_minus_theta = ROOT.TH1F("kaon_minus_theta", "Kaon- Theta", 500, 0, 50)
    hist_W_values = ROOT.TH1F("W_values", "W", 500, 0, 5)
    hist_Q2_values = ROOT.TH1F("Q2_values", "Q²", 500, 0, 15)
    hist_x_values = ROOT.TH1F("x_values", "x", 500, 0, 1)
    hist_t_values = ROOT.TH1F("t_values", "t", 500, 0, 20)
    hist_filtered_W_values = ROOT.TH1F("filtered_W_values", "Reconstructed W", 500, 0, 5)
    hist_filtered_Q2_values = ROOT.TH1F("filtered_Q2_values", "Reconstructed Q²", 500, 0, 15)
    hist_filtered_x_values = ROOT.TH1F("filtered_x_values", "Reconstructed x", 500, 0, 1)
    hist_filtered_t_values = ROOT.TH1F("filtered_t_values", "Reconstructed t", 500, 0, 20)
    hist_Q2_resolutions = ROOT.TH1F("Q2_resolutions", "Q² Resolutions", 500, -5, 5)
    hist_W_resolutions = ROOT.TH1F("W_resolutions", "W Resolutions", 500, -2, 2)
    hist_x_resolutions = ROOT.TH1F("x_resolutions", "x Resolutions", 500, -1, 1)
    hist_y_resolutions = ROOT.TH1F("y_resolutions", "y Resolutions", 500, -1, 1)
    hist_t_resolutions = ROOT.TH1F("t_resolutions", "t Resolutions", 500, -5, 5)
    hist_missing_mass = ROOT.TH1F("missing_mass", "Missing Mass", 500, -1, 2)
    hist_KK_mass = ROOT.TH1F("KK_mass", "K+ K- Invariant Mass", 500, 0, 3)

    for value in electron_momenta:
        hist_electron_momenta.Fill(value)
    for value in proton_momenta:
        hist_proton_momenta.Fill(value)
    for value in kaon_plus_momenta:
        hist_kaon_plus_momenta.Fill(value)
    for value in kaon_minus_momenta:
        hist_kaon_minus_momenta.Fill(value)
    for value in electron_theta:
        hist_electron_theta.Fill(value)
    for value in proton_theta:
        hist_proton_theta.Fill(value)
    for value in kaon_plus_theta:
        hist_kaon_plus_theta.Fill(value)
    for value in kaon_minus_theta:
        hist_kaon_minus_theta.Fill(value)
    for value in W_values:
        hist_W_values.Fill(value)
    for value in Q2_values:
        hist_Q2_values.Fill(value)
    for value in x_values:
        hist_x_values.Fill(value)
    for value in t_values:
        hist_t_values.Fill(value)
    for value in filtered_W_values:
        hist_filtered_W_values.Fill(value)
    for value in filtered_Q2_values:
        hist_filtered_Q2_values.Fill(value)
    for value in filtered_x_values:
        hist_filtered_x_values.Fill(value)
    for value in filtered_t_values:
        hist_filtered_t_values.Fill(value)
    for value in Q2_resolutions:
        hist_Q2_resolutions.Fill(value)
    for value in W_resolutions:
        hist_W_resolutions.Fill(value)
    for value in x_resolutions:
        hist_x_resolutions.Fill(value)
    for value in y_resolutions:
        hist_y_resolutions.Fill(value)
    for value in t_resolutions:
        hist_t_resolutions.Fill(value)
    for value in missing_mass_values:
        hist_missing_mass.Fill(value)
    for value in KK_mass_values:
        hist_KK_mass.Fill(value)

    hist_electron_momenta.Write()
    hist_proton_momenta.Write()
    hist_kaon_plus_momenta.Write()
    hist_kaon_minus_momenta.Write()
    hist_electron_theta.Write()
    hist_proton_theta.Write()
    hist_kaon_plus_theta.Write()
    hist_kaon_minus_theta.Write()
    hist_W_values.Write()
    hist_Q2_values.Write()
    hist_x_values.Write()
    hist_t_values.Write()
    hist_filtered_W_values.Write()
    hist_filtered_Q2_values.Write()
    hist_filtered_x_values.Write()
    hist_filtered_t_values.Write()
    hist_Q2_resolutions.Write()
    hist_W_resolutions.Write()
    hist_x_resolutions.Write()
    hist_y_resolutions.Write()
    hist_t_resolutions.Write()
    hist_missing_mass.Write()
    hist_KK_mass.Write()

    output_root_file.Close()

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
