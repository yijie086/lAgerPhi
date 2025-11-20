import ROOT
import numpy as np
import sys
import random

def initialize_histograms(hist_prefix):
    histograms = {
        f"{hist_prefix}_electron_momenta": ROOT.TH1F(f"{hist_prefix}_electron_momenta", "Electron Momenta", 500, 0, 10.6),
        f"{hist_prefix}_proton_momenta": ROOT.TH1F(f"{hist_prefix}_proton_momenta", "Proton Momenta", 500, 0, 10.6),
        f"{hist_prefix}_kaon_plus_momenta": ROOT.TH1F(f"{hist_prefix}_kaon_plus_momenta", "Kaon+ Momenta", 500, 0, 10.6),
        f"{hist_prefix}_kaon_minus_momenta": ROOT.TH1F(f"{hist_prefix}_kaon_minus_momenta", "Kaon- Momenta", 500, 0, 10.6),
        f"{hist_prefix}_electron_theta": ROOT.TH1F(f"{hist_prefix}_electron_theta", "Electron Theta", 500, 0, 50),
        f"{hist_prefix}_proton_theta": ROOT.TH1F(f"{hist_prefix}_proton_theta", "Proton Theta", 500, 0, 50),
        f"{hist_prefix}_kaon_plus_theta": ROOT.TH1F(f"{hist_prefix}_kaon_plus_theta", "Kaon+ Theta", 500, 0, 50),
        f"{hist_prefix}_kaon_minus_theta": ROOT.TH1F(f"{hist_prefix}_kaon_minus_theta", "Kaon- Theta", 500, 0, 50),
        f"{hist_prefix}_W_values": ROOT.TH1F(f"{hist_prefix}_W_values", "W", 500, 0, 5),
        f"{hist_prefix}_Q2_values": ROOT.TH1F(f"{hist_prefix}_Q2_values", "Q²", 500, 0, 15),
        f"{hist_prefix}_x_values": ROOT.TH1F(f"{hist_prefix}_x_values", "x", 500, 0, 1),
        f"{hist_prefix}_t_values": ROOT.TH1F(f"{hist_prefix}_t_values", "t", 500, 0, 20),
        f"{hist_prefix}_filtered_W_values": ROOT.TH1F(f"{hist_prefix}_filtered_W_values", "Reconstructed W", 500, 0, 5),
        f"{hist_prefix}_filtered_Q2_values": ROOT.TH1F(f"{hist_prefix}_filtered_Q2_values", "Reconstructed Q²", 500, 0, 15),
        f"{hist_prefix}_filtered_x_values": ROOT.TH1F(f"{hist_prefix}_filtered_x_values", "Reconstructed x", 500, 0, 1),
        f"{hist_prefix}_filtered_t_values": ROOT.TH1F(f"{hist_prefix}_filtered_t_values", "Reconstructed t", 500, 0, 20),
        f"{hist_prefix}_Q2_resolutions": ROOT.TH1F(f"{hist_prefix}_Q2_resolutions", "Q² Resolutions", 500, -5, 5),
        f"{hist_prefix}_W_resolutions": ROOT.TH1F(f"{hist_prefix}_W_resolutions", "W Resolutions", 500, -2, 2),
        f"{hist_prefix}_x_resolutions": ROOT.TH1F(f"{hist_prefix}_x_resolutions", "x Resolutions", 500, -1, 1),
        f"{hist_prefix}_y_resolutions": ROOT.TH1F(f"{hist_prefix}_y_resolutions", "y Resolutions", 500, -1, 1),
        f"{hist_prefix}_t_resolutions": ROOT.TH1F(f"{hist_prefix}_t_resolutions", "t Resolutions", 500, -5, 5),
        f"{hist_prefix}_missing_mass": ROOT.TH1F(f"{hist_prefix}_missing_mass", "Missing Mass", 500, -1, 2),
        f"{hist_prefix}_KK_mass": ROOT.TH1F(f"{hist_prefix}_KK_mass", "K+ K- Invariant Mass", 500, 0, 3),
    }
    return histograms

def process_events(input_file, hist_prefix, histograms):
    # Open the ROOT file
    file = ROOT.TFile(input_file)
    tree = None
    # Get the TTree
    if hist_prefix == "sig":  # Running over Signal
        tree = file.Get("lAger")
    if hist_prefix == "bk":  # Running over background
        tree = file.Get("DiffPhiBkgdEvents")
    
    if not tree:
        print(f"Error: TTree not found in {input_file}")
        return {}

    print(f"Processing TTree {tree.GetName()} with {tree.GetEntries()} entries in {input_file}")

    # Cuts (modify these values as needed)
    W_cut = (1.0, 2.4)
    Q2_cut = (6.0, 60.0)
    x_cut = (0.0001, 0.9999)

    # Define initial state four-momentum (proton at rest, incoming electron)
    initial_proton = ROOT.TLorentzVector(0, 0, 0, 0.938)  # Proton mass in GeV
    initial_electron = ROOT.TLorentzVector(0, 0, 10.6, np.sqrt(10.6**2 + 0.000511**2))  # Incoming electron 10.6 GeV in z-direction

    # Function to calculate polar angle theta in degrees
    def calculate_theta(p):
        if p.Pz() == 0:
            return 0
        else:
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

    def create_lorentz_vector(particle):
        return ROOT.TLorentzVector(particle.Px(), particle.Py(), particle.Pz(), particle.Energy())

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
    electron = None
    proton = None
    kaon_plus = None
    kaon_minus = None

    # Loop over the entries in the tree
    event_no = 0
    for event in tree:
        if event_no%100000 == 0:
            print(f"Processing Event Number {event_no}")
        
        event_no = event_no+1
        
        if hist_prefix == "bk":
            electron = event.Electron
            charged_particles = [
                (event.ChargedParticle1, event.ChargedPID1),
                (event.ChargedParticle2, event.ChargedPID2),
                (event.ChargedParticle3, event.ChargedPID3),
                (event.ChargedParticle4, event.ChargedPID4),
                (event.ChargedParticle5, event.ChargedPID5),
            ]
           
            for charged_particle, pid in charged_particles:
                if pid == 2212:  # Proton
                    proton = charged_particle
                elif pid == -321:  # K+
                    kaon_plus = charged_particle
                elif pid == 321:  # K-
                    kaon_minus = charged_particle
                elif pid == 211:
                    if random.random() < 0.02:  # 2% chance
                        kaon_plus = charged_particle
                elif pid == -211:
                    if random.random() < 0.02:  # 2% chance
                        kaon_minus = charged_particle
                        
            if not proton or not kaon_plus or not kaon_minus:
                continue
        elif hist_prefix == "sig":
            particles = event.particles
            electron = create_lorentz_vector(particles[4])
            proton = create_lorentz_vector(particles[6])
            kaon_minus = create_lorentz_vector(particles[7])
            kaon_plus = create_lorentz_vector(particles[8])
            
        # Calculate kinematic variables
        W, Q2, x, y = calculate_kinematics(electron)
        t = (initial_proton - proton).M2()
        #print(f"x = {x}, Q2 = {Q2}, W = {W}, t = {t} ")
        # Apply generator-level cuts
        if (W_cut[0] <= W <= W_cut[1]) and (Q2_cut[0] <= Q2 <= Q2_cut[1]) and (x_cut[0] <= x <= x_cut[1]):
            # Append generator-level kinematic variables
            #print(f"inside generator cuts x = {x}, Q2 = {Q2}, W = {W}, t = {t} ")
            histograms[f"{hist_prefix}_W_values"].Fill(W)
            histograms[f"{hist_prefix}_Q2_values"].Fill(Q2)
            histograms[f"{hist_prefix}_x_values"].Fill(x)
            histograms[f"{hist_prefix}_t_values"].Fill(t)

            # Append generator-level momenta and polar angles for all particles
            histograms[f"{hist_prefix}_electron_momenta"].Fill(electron.P())
            histograms[f"{hist_prefix}_proton_momenta"].Fill(proton.P())
            histograms[f"{hist_prefix}_kaon_plus_momenta"].Fill(kaon_plus.P())
            histograms[f"{hist_prefix}_kaon_minus_momenta"].Fill(kaon_minus.P())

            histograms[f"{hist_prefix}_electron_theta"].Fill(calculate_theta(electron))
            histograms[f"{hist_prefix}_proton_theta"].Fill(calculate_theta(proton))
            histograms[f"{hist_prefix}_kaon_plus_theta"].Fill(calculate_theta(kaon_plus))
            histograms[f"{hist_prefix}_kaon_minus_theta"].Fill(calculate_theta(kaon_minus))

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
        proton_bin_x = overall_acceptance_hist.GetXaxis().FindBin(proton_theta_val)
        proton_bin_y = overall_acceptance_hist.GetYaxis().FindBin(proton_momentum_val)
        proton_acceptance = overall_acceptance_hist.GetBinContent(proton_bin_x, proton_bin_y)
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

            filtered_W, filtered_Q2, filtered_x, filtered_y = calculate_kinematics(smeared_electron)

            smeared_proton_momentum = smear_momentum(proton.P())
            smeared_proton_theta = smear_angle(np.radians(calculate_theta(proton)))
            smeared_proton_phi = smear_angle(proton.Phi())

            smeared_proton_px = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.cos(smeared_proton_phi)
            smeared_proton_py = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.sin(smeared_proton_phi)
            smeared_proton_pz = smeared_proton_momentum * np.cos(smeared_proton_theta)
            smeared_proton_E = np.sqrt(smeared_proton_px**2 + smeared_proton_py**2 + smeared_proton_pz**2 + 0.938**2)
            smeared_proton = ROOT.TLorentzVector()
            smeared_proton.SetPxPyPzE(smeared_proton_px, smeared_proton_py, smeared_proton_pz, smeared_proton_E)

            filtered_t = (proton - initial_proton).M2()

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

            if (W_cut[0] <= filtered_W <= W_cut[1]) and (Q2_cut[0] <= filtered_Q2 <= Q2_cut[1]) and (x_cut[0] <= filtered_x <= x_cut[1]):
                #print(f"Reco Event {event_no} passed acceptance cuts")
                histograms[f"{hist_prefix}_filtered_W_values"].Fill(filtered_W)
                histograms[f"{hist_prefix}_filtered_Q2_values"].Fill(filtered_Q2)
                histograms[f"{hist_prefix}_filtered_x_values"].Fill(filtered_x)
                histograms[f"{hist_prefix}_filtered_t_values"].Fill(filtered_t)
                
                # Calculate resolutions
                Q2_resolution = filtered_Q2 - Q2
                W_resolution = filtered_W - W
                x_resolution = filtered_x - x
                y_resolution = filtered_y - y
                t_resolution = (initial_proton - smeared_proton).M2() - t
                
                histograms[f"{hist_prefix}_Q2_resolutions"].Fill(Q2_resolution)
                histograms[f"{hist_prefix}_W_resolutions"].Fill(W_resolution)
                histograms[f"{hist_prefix}_x_resolutions"].Fill(x_resolution)
                histograms[f"{hist_prefix}_y_resolutions"].Fill(y_resolution)
                histograms[f"{hist_prefix}_t_resolutions"].Fill(t_resolution)
                reconstructed_particles = [smeared_electron, smeared_kaon_minus, smeared_kaon_plus, smeared_proton]
                
                missing_mass = ((initial_electron + initial_proton) - sum(reconstructed_particles, ROOT.TLorentzVector())).M()
                histograms[f"{hist_prefix}_missing_mass"].Fill(missing_mass)
                
                KK_mass = (smeared_kaon_plus+smeared_kaon_minus).M()
                histograms[f"{hist_prefix}_KK_mass"].Fill(KK_mass)
                
    return histograms

def main(signal_file, output_file, background_file=None):
    # Initialize histograms
    signal_histograms = initialize_histograms("sig")
    background_histograms = initialize_histograms("bk") if background_file else {}

    # Process signal file
    process_events(signal_file, "sig", signal_histograms)

    # Process background file if provided
    if background_file:
        process_events(background_file, "bk", background_histograms)

    # Save histograms to ROOT file
    output_root_file = ROOT.TFile(output_file, "RECREATE")

    # Write signal histograms
    print("Writing signal histograms...")
    for hist_name, hist in signal_histograms.items():
        if hist:  # Ensure the histogram is not None
            hist.Write()
        else:
            print(f"Skipping {hist_name} because it is None")

    # Write background histograms if provided
    if background_file:
        print("Writing background histograms...")
        for hist_name, hist in background_histograms.items():
            if hist:  # Ensure the histogram is not None
                hist.Write()
            else:
                print(f"Skipping {hist_name} because it is None")

    output_root_file.Close()
    print(f"Output file {output_file} closed")
    
if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        print("Usage: python script.py <signal_file> <output_file> [<background_file>]")
        sys.exit(1)

    signal_file = sys.argv[1]
    output_file = sys.argv[2]
    background_file = sys.argv[3] if len(sys.argv) == 4 else None

    main(signal_file, output_file, background_file)
