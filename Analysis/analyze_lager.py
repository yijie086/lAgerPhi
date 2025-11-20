import ROOT
import numpy as np
import matplotlib.pyplot as plt
import time

# Define the file paths
file_path = "/work/hallc/jpsi-007/klest/lager2/outputs/GenOnly-ep-phi.ep-phi.4pi.run01001-lumi100.root"
acceptance_file_path = "/home/klest/pcsim_proc_detector_solid-fastmc_acceptance_solid_JPsi_electron_target315_output.root"

# Open the ROOT file
file = ROOT.TFile(file_path)

# Get the TTree
tree = file.Get("lAger")

# Initialize lists to store momenta and polar angles
electron_momenta = []
phi_momenta = []
proton_momenta = []
kaon1_momenta = []
kaon2_momenta = []

electron_theta = []
phi_theta = []
proton_theta = []
kaon1_theta = []
kaon2_theta = []

# Kinematic variables
W_values = []
Q2_values = []
x_values = []
t_values = []

filtered_electron_momenta = []
filtered_proton_momenta = []
filtered_kaon1_momenta = []
filtered_kaon2_momenta = []

filtered_electron_theta = []
filtered_proton_theta = []
filtered_kaon1_theta = []
filtered_kaon2_theta = []

# Kinematic variables after filtering
filtered_W_values = []
filtered_Q2_values = []
filtered_x_values = []
filtered_t_values = []

# Missing mass values for different cases
missing_mass_kaon_missing = []
missing_mass_proton_missing = []

# Resolutions
Q2_resolutions = []
W_resolutions = []
x_resolutions = []
t_resolutions = []
y_resolutions = []

# Cuts (modify these values as needed)
W_cut = (1.0, 2.4)
Q2_cut = (5, 60)
x_cut = (0.0001, 0.9999)

# Define initial state four-momentum (proton at rest, incoming electron)
initial_proton = ROOT.TLorentzVector(0, 0, 0, 0.938)  # Proton mass in GeV
initial_electron = ROOT.TLorentzVector(0, 0, 10.6, np.sqrt(10.6**2 + 0.000511**2))  # Incoming electron 10.6 GeV in z-direction

# Function to calculate polar angle theta in degrees
def calculate_theta(p):
    return np.degrees(np.arccos(p.Pz() / p.P()))

# Function to create a TLorentzVector from a TParticle
def create_lorentz_vector(particle):
    return ROOT.TLorentzVector(particle.Px(), particle.Py(), particle.Pz(), particle.Energy())

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
acceptance_file = ROOT.TFile(acceptance_file_path)
kaon_acceptance_hist = acceptance_file.Get("acceptance_ThetaP_forwardangle")
overall_acceptance_hist = acceptance_file.Get("acceptance_ThetaP_overall")

def scale_histogram(hist, max_value):
    hist_max = hist.GetMaximum()
    scale_factor = max_value / hist_max
    hist.Scale(scale_factor)

# Scale the acceptance histograms
scale_histogram(kaon_acceptance_hist, 0.9)
scale_histogram(overall_acceptance_hist, 0.9)
start_time = time.time()
iteration_count = 0
full_event_meas_count = 0
most_event_meas_count = 0

# Loop over the entries in the tree
for event in tree:
    particles = event.particles
    if iteration_count % 100000 == 0:
        print(f"Processing event number: {iteration_count}")

    iteration_count += 1
    # Get the relevant particles
    electron = create_lorentz_vector(particles[4])
    proton = create_lorentz_vector(particles[6])
    kaon1 = create_lorentz_vector(particles[7])
    kaon2 = create_lorentz_vector(particles[8])

    # Calculate kinematic variables
    W, Q2, x, y = calculate_kinematics(electron)
    t = (proton - initial_proton).M2()
    # Apply cuts
    if not (W_cut[0] <= W <= W_cut[1]):
        continue
    if not (Q2_cut[0] <= Q2 <= Q2_cut[1]):
        continue
    if not (x_cut[0] <= x <= x_cut[1]):
        continue

    # Append kinematic variables
    W_values.append(W)
    Q2_values.append(Q2)
    x_values.append(x)
    t_values.append(t)

    # Append momenta and polar angles for all particles
    electron_momenta.append(electron.P())
    proton_momenta.append(proton.P())
    kaon1_momenta.append(kaon1.P())
    kaon2_momenta.append(kaon2.P())

    electron_theta.append(calculate_theta(electron))
    proton_theta.append(calculate_theta(proton))
    kaon1_theta.append(calculate_theta(kaon1))
    kaon2_theta.append(calculate_theta(kaon2))

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

    # Kaon1 acceptance
    kaon1_theta_val = calculate_theta(kaon1)
    kaon1_momentum_val = kaon1.P()
    kaon1_bin_x = kaon_acceptance_hist.GetXaxis().FindBin(kaon1_theta_val)
    kaon1_bin_y = kaon_acceptance_hist.GetYaxis().FindBin(kaon1_momentum_val)
    kaon1_acceptance = kaon_acceptance_hist.GetBinContent(kaon1_bin_x, kaon1_bin_y)
    if np.random.random() < kaon1_acceptance:
        acceptances.append(True)
    else:
        acceptances.append(False)

    # Kaon2 acceptance
    kaon2_theta_val = calculate_theta(kaon2)
    kaon2_momentum_val = kaon2.P()
    kaon2_bin_x = kaon_acceptance_hist.GetXaxis().FindBin(kaon2_theta_val)
    kaon2_bin_y = kaon_acceptance_hist.GetYaxis().FindBin(kaon2_momentum_val)
    kaon2_acceptance = kaon_acceptance_hist.GetBinContent(kaon2_bin_x, kaon2_bin_y)
    if np.random.random() < kaon2_acceptance:
        acceptances.append(True)
    else:
        acceptances.append(False)

    # Event passes if electron and at least two of the three (kaon1, kaon2, proton) pass
    if acceptances[0] and sum(acceptances[1:]) >= 2:
        most_event_meas_count += 1
        smeared_electron_momentum = smear_momentum(electron.P())
        smeared_electron_theta = smear_angle(np.radians(calculate_theta(electron)))
        smeared_electron_phi = smear_angle(electron.Phi())

        smeared_electron_px = smeared_electron_momentum * np.sin(smeared_electron_theta) * np.cos(smeared_electron_phi)
        smeared_electron_py = smeared_electron_momentum * np.sin(smeared_electron_theta) * np.sin(smeared_electron_phi)
        smeared_electron_pz = smeared_electron_momentum * np.cos(smeared_electron_theta)

        smeared_electron = ROOT.TLorentzVector()
        smeared_electron.SetPxPyPzE(smeared_electron_px, smeared_electron_py, smeared_electron_pz, electron.E())

        filtered_electron_momenta.append(smeared_electron.P())
        filtered_electron_theta.append(np.degrees(smeared_electron_theta))
        filtered_W, filtered_Q2, filtered_x, filtered_y = calculate_kinematics(smeared_electron)

        filtered_W_values.append(filtered_W)
        filtered_Q2_values.append(filtered_Q2)
        filtered_x_values.append(filtered_x)

        # Calculate resolutions
        Q2_resolution = filtered_Q2 - Q2
        W_resolution = filtered_W - W
        x_resolution = filtered_x - x
        y_resolution = filtered_y - y
        

        Q2_resolutions.append(Q2_resolution)
        W_resolutions.append(W_resolution)
        x_resolutions.append(x_resolution)
        y_resolutions.append(y_resolution)
        

        if acceptances[1]:  # Proton
            smeared_proton_momentum = smear_momentum(proton.P())
            smeared_proton_theta = smear_angle(np.radians(calculate_theta(proton)))
            smeared_proton_phi = smear_angle(proton.Phi())

            smeared_proton_px = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.cos(smeared_proton_phi)
            smeared_proton_py = smeared_proton_momentum * np.sin(smeared_proton_theta) * np.sin(smeared_proton_phi)
            smeared_proton_pz = smeared_proton_momentum * np.cos(smeared_proton_theta)

            smeared_proton = ROOT.TLorentzVector()
            smeared_proton.SetPxPyPzE(smeared_proton_px, smeared_proton_py, smeared_proton_pz, proton.E())
            t_resolution = (smeared_proton - initial_proton).M2() - t
            t_resolutions.append(t_resolution)
            filtered_proton_momenta.append(smeared_proton.P())
            filtered_proton_theta.append(np.degrees(smeared_proton_theta))
            filtered_t = (smeared_proton - initial_proton).M2()
            filtered_t_values.append(filtered_t)

        if acceptances[2]:  # Kaon1
            smeared_kaon1_momentum = smear_momentum(kaon1.P())
            smeared_kaon1_theta = smear_angle(np.radians(calculate_theta(kaon1)))
            smeared_kaon1_phi = smear_angle(kaon1.Phi())

            smeared_kaon1_px = smeared_kaon1_momentum * np.sin(smeared_kaon1_theta) * np.cos(smeared_kaon1_phi)
            smeared_kaon1_py = smeared_kaon1_momentum * np.sin(smeared_kaon1_theta) * np.sin(smeared_kaon1_phi)
            smeared_kaon1_pz = smeared_kaon1_momentum * np.cos(smeared_kaon1_theta)

            smeared_kaon1 = ROOT.TLorentzVector()
            smeared_kaon1.SetPxPyPzE(smeared_kaon1_px, smeared_kaon1_py, smeared_kaon1_pz, kaon1.E())

            filtered_kaon1_momenta.append(smeared_kaon1.P())
            filtered_kaon1_theta.append(np.degrees(smeared_kaon1_theta))

        if acceptances[3]:  # Kaon2
            smeared_kaon2_momentum = smear_momentum(kaon2.P())
            smeared_kaon2_theta = smear_angle(np.radians(calculate_theta(kaon2)))
            smeared_kaon2_phi = smear_angle(kaon2.Phi())

            smeared_kaon2_px = smeared_kaon2_momentum * np.sin(smeared_kaon2_theta) * np.cos(smeared_kaon2_phi)
            smeared_kaon2_py = smeared_kaon2_momentum * np.sin(smeared_kaon2_theta) * np.sin(smeared_kaon2_phi)
            smeared_kaon2_pz = smeared_kaon2_momentum * np.cos(smeared_kaon2_theta)

            smeared_kaon2 = ROOT.TLorentzVector()
            smeared_kaon2.SetPxPyPzE(smeared_kaon2_px, smeared_kaon2_py, smeared_kaon2_pz, kaon2.E())

            filtered_kaon2_momenta.append(smeared_kaon2.P())
            filtered_kaon2_theta.append(np.degrees(smeared_kaon2_theta))

        # Calculate missing mass with smeared particles
        if not acceptances[1]:  # Proton missed
            reconstructed_particles = [smeared_electron, smeared_kaon1, smeared_kaon2]
            missing_mass = ((initial_electron + initial_proton) - sum(reconstructed_particles, ROOT.TLorentzVector())).M()
            missing_mass_proton_missing.append(missing_mass)
        elif not (acceptances[2] and acceptances[3]):  # One of the Kaons missed
            reconstructed_particles = [smeared_electron, smeared_proton]
            if acceptances[2]:
                reconstructed_particles.append(smeared_kaon1)
            if acceptances[3]:
                reconstructed_particles.append(smeared_kaon2)
            missing_mass = ((initial_electron + initial_proton) - sum(reconstructed_particles, ROOT.TLorentzVector())).M()
            missing_mass_kaon_missing.append(missing_mass)

        if acceptances[1] and acceptances[2] and acceptances[3]:  # all reconstructed
            full_event_meas_count += 1

total_elapsed_time = time.time() - start_time
print(f"Event loop completed in {total_elapsed_time:.2f} seconds!")
print(f"Total number of events is {iteration_count}, number of fully reconstructed events is {full_event_meas_count}, full acceptance ratio is {full_event_meas_count/iteration_count}, number of 2/3 reconstructed is {most_event_meas_count}, ratio is {most_event_meas_count/iteration_count}")

# Define binning parameters
momentum_bins = np.linspace(0, 10.6, 50)  # JLab range for momenta
theta_bins = np.linspace(0, 50, 50)       # JLab range for polar angle theta
mx_bins = np.linspace(-1, 1.5, 50)        # JLab range for polar angle theta
resolution_bins = np.linspace(-0.5, 0.5, 50)

# Create a ROOT file to save histograms
output_root_file = ROOT.TFile("output_histograms.root", "RECREATE")

# Plotting
plt.figure(figsize=(18, 18))

# Momentum distribution
plt.hist(electron_momenta, bins=momentum_bins, alpha=0.5, label='Electron', color='blue')
plt.hist(proton_momenta, bins=momentum_bins, alpha=0.5, label='Proton', color='green')
plt.hist(kaon1_momenta, bins=momentum_bins, alpha=0.5, label='$K^-$', color='red')
plt.hist(kaon2_momenta, bins=momentum_bins, alpha=0.5, label='$K^+$', color='orange')
plt.xlabel('Momentum (GeV)')
plt.ylabel('Counts')
plt.legend()
plt.title('Momentum Distribution')
plt.savefig('momentum_distribution.pdf')

# Theta distribution
plt.hist(electron_theta, bins=theta_bins, alpha=0.5, label='Electron', color='blue')
plt.hist(proton_theta, bins=theta_bins, alpha=0.5, label='Proton', color='green')
plt.hist(kaon1_theta, bins=theta_bins, alpha=0.5, label='$K^-$', color='red')
plt.hist(kaon2_theta, bins=theta_bins, alpha=0.5, label='$K^+$', color='orange')
plt.xlabel('Theta (degrees)')
plt.ylabel('Counts')
plt.legend()
plt.title('Theta Distribution')
plt.savefig('theta_distribution.pdf')

plt.figure(figsize=(18, 12))

# W Distribution
plt.subplot(2, 2, 1)
plt.hist(W_values, bins=50, alpha=0.5)
plt.xlabel('W (GeV)')
plt.ylabel('Counts')
plt.title('W Distribution')

# Q² Distribution
plt.subplot(2, 2, 2)
plt.hist(Q2_values, bins=50, alpha=0.5)
plt.xlabel('Q² (GeV²)')
plt.ylabel('Counts')
plt.title('Q² Distribution')

# x Distribution
plt.subplot(2, 2, 3)
plt.hist(x_values, bins=50, alpha=0.5)
plt.xlabel('x')
plt.ylabel('Counts')
plt.title('x Distribution')

# t Distribution
plt.subplot(2, 2, 4)
plt.hist(t_values, bins=50, alpha=0.5)
plt.xlabel('t (GeV²)')
plt.ylabel('Counts')
plt.title('t Distribution')

plt.tight_layout()
plt.savefig('kinematic_variables_distribution.pdf')

plt.figure(figsize=(18, 12))

# Filtered Theta vs Momentum (Kaon1)
plt.subplot(2, 2, 1)
plt.hist2d(filtered_kaon1_theta, filtered_kaon1_momenta, bins=[theta_bins, momentum_bins], cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('Theta (degrees)')
plt.ylabel('Momentum (GeV)')
plt.title('$K^-$ Theta vs Momentum (Filtered by Acceptance)')

# Filtered Theta vs Momentum (Electron)
plt.subplot(2, 2, 2)
plt.hist2d(filtered_electron_theta, filtered_electron_momenta, bins=[theta_bins, momentum_bins], cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('Theta (degrees)')
plt.ylabel('Momentum (GeV)')
plt.title('Electron Theta vs Momentum (Filtered by Acceptance)')

# Missing Mass Distribution
plt.subplot(2, 2, 3)
plt.hist(missing_mass_kaon_missing, bins=mx_bins, alpha=0.5, color='red', label='Kaon Missing')
plt.hist(missing_mass_proton_missing, bins=mx_bins, alpha=0.5, color='green', label='Proton Missing')
plt.xlabel('Missing Mass (GeV)')
plt.ylabel('Counts')
plt.title('Missing Mass Distribution')
plt.legend()

plt.tight_layout()
plt.savefig('missing_mass_distribution.pdf')

# Separate canvas for kinematic variables after acceptance
plt.figure(figsize=(24, 6))

plt.subplot(1, 4, 1)
plt.hist(filtered_W_values, bins=50, alpha=0.5)
plt.xlabel('W (GeV)')
plt.ylabel('Counts')
plt.title('W Distribution (After Acceptance)')

plt.subplot(1, 4, 2)
plt.hist(filtered_Q2_values, bins=50, alpha=0.5)
plt.xlabel('Q² (GeV²)')
plt.ylabel('Counts')
plt.title('Q² Distribution (After Acceptance)')

plt.subplot(1, 4, 3)
plt.hist(filtered_x_values, bins=50, alpha=0.5)
plt.xlabel('x')
plt.ylabel('Counts')
plt.title('x Distribution (After Acceptance)')

plt.subplot(1, 4, 4)
plt.hist(filtered_t_values, bins=50, alpha=0.5)
plt.xlabel('t (GeV²)')
plt.ylabel('Counts')
plt.title('t Distribution (After Acceptance)')

plt.tight_layout()
plt.savefig('filtered_kinematic_variables_distribution.pdf')

# Plotting resolutions
plt.figure(figsize=(18, 12))

plt.subplot(2, 2, 1)
plt.hist(Q2_resolutions, bins=resolution_bins, alpha=0.5, label='Q² resolution', color='blue')
plt.xlabel('Reco Q² - Gen Q²')
plt.ylabel('Counts')
plt.title('Q² Resolution')
plt.legend()

plt.subplot(2, 2, 2)
plt.hist(W_resolutions, bins=resolution_bins, alpha=0.5, label='W resolution', color='green')
plt.xlabel('Reco W - Gen W')
plt.ylabel('Counts')
plt.title('W Resolution')
plt.legend()

plt.subplot(2, 2, 3)
plt.hist(x_resolutions, bins=resolution_bins, alpha=0.5, label='x resolution', color='red')
plt.xlabel('Reco x - Gen x')
plt.ylabel('Counts')
plt.title('x Resolution')
plt.legend()

plt.subplot(2, 2, 4)
plt.hist(t_resolutions, bins=resolution_bins, alpha=0.5, label='t resolution', color='orange')
plt.xlabel('Reco t - Gen t')
plt.ylabel('Counts')
plt.title('t Resolution')
plt.legend()

plt.tight_layout()
plt.savefig('resolutions.pdf')

# Save histograms to ROOT file
output_root_file.cd()

hist_electron_momenta = ROOT.TH1F("electron_momenta", "Electron Momenta", len(momentum_bins)-1, momentum_bins)
hist_proton_momenta = ROOT.TH1F("proton_momenta", "Proton Momenta", len(momentum_bins)-1, momentum_bins)
hist_kaon1_momenta = ROOT.TH1F("kaon1_momenta", "K^- Momenta", len(momentum_bins)-1, momentum_bins)
hist_kaon2_momenta = ROOT.TH1F("kaon2_momenta", "K^+ Momenta", len(momentum_bins)-1, momentum_bins)
hist_electron_theta = ROOT.TH1F("electron_theta", "Electron Theta", len(theta_bins)-1, theta_bins)
hist_proton_theta = ROOT.TH1F("proton_theta", "Proton Theta", len(theta_bins)-1, theta_bins)
hist_kaon1_theta = ROOT.TH1F("kaon1_theta", "K^- Theta", len(theta_bins)-1, theta_bins)
hist_kaon2_theta = ROOT.TH1F("kaon2_theta", "K^+ Theta", len(theta_bins)-1, theta_bins)
hist_W_values = ROOT.TH1F("W_values", "W", 50, min(W_values), max(W_values))
hist_Q2_values = ROOT.TH1F("Q2_values", "Q²", 50, min(Q2_values), max(Q2_values))
hist_x_values = ROOT.TH1F("x_values", "x", 50, min(x_values), max(x_values))
hist_t_values = ROOT.TH1F("t_values", "t", 50, min(t_values), max(t_values))
hist_Q2_resolutions = ROOT.TH1F("Q2_resolutions", "Q² Resolutions", len(resolution_bins)-1, resolution_bins)
hist_W_resolutions = ROOT.TH1F("W_resolutions", "W Resolutions", len(resolution_bins)-1, resolution_bins)
hist_x_resolutions = ROOT.TH1F("x_resolutions", "x Resolutions", len(resolution_bins)-1, resolution_bins)
hist_t_resolutions = ROOT.TH1F("t_resolutions", "t Resolutions", len(resolution_bins)-1, resolution_bins)

for value in electron_momenta:
    hist_electron_momenta.Fill(value)
for value in proton_momenta:
    hist_proton_momenta.Fill(value)
for value in kaon1_momenta:
    hist_kaon1_momenta.Fill(value)
for value in kaon2_momenta:
    hist_kaon2_momenta.Fill(value)
for value in electron_theta:
    hist_electron_theta.Fill(value)
for value in proton_theta:
    hist_proton_theta.Fill(value)
for value in kaon1_theta:
    hist_kaon1_theta.Fill(value)
for value in kaon2_theta:
    hist_kaon2_theta.Fill(value)
for value in W_values:
    hist_W_values.Fill(value)
for value in Q2_values:
    hist_Q2_values.Fill(value)
for value in x_values:
    hist_x_values.Fill(value)
for value in t_values:
    hist_t_values.Fill(value)
for value in Q2_resolutions:
    hist_Q2_resolutions.Fill(value)
for value in W_resolutions:
    hist_W_resolutions.Fill(value)
for value in x_resolutions:
    hist_x_resolutions.Fill(value)
for value in t_resolutions:
    hist_t_resolutions.Fill(value)

hist_electron_momenta.Write()
hist_proton_momenta.Write()
hist_kaon1_momenta.Write()
hist_kaon2_momenta.Write()
hist_electron_theta.Write()
hist_proton_theta.Write()
hist_kaon1_theta.Write()
hist_kaon2_theta.Write()
hist_W_values.Write()
hist_Q2_values.Write()
hist_x_values.Write()
hist_t_values.Write()
hist_Q2_resolutions.Write()
hist_W_resolutions.Write()
hist_x_resolutions.Write()
hist_t_resolutions.Write()

output_root_file.Close()
