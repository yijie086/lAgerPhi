import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile("outputs/EIC-ep-phi.ep-phi.4pi.run00888-lumi5.root")

# Get the TTree
tree = file.Get("lAger")

# Initialize lists to store momenta and pseudorapidity
electron_momenta = []
phi_momenta = []
proton_momenta = []
kaon1_momenta = []
kaon2_momenta = []

electron_eta = []
phi_eta = []
proton_eta = []
kaon1_eta = []
kaon2_eta = []

# Kinematic variables
W_values = []
Q2_values = []
x_values = []

# Cuts (modify these values as needed)
W_cut = (1.96, 3.0)
Q2_cut = (1.0, 100000.0)
x_cut = (0.00001, 0.9999)

# Define initial state four-momentum (proton at 41 GeV in +z direction, electron at 5 GeV in -z direction)
initial_proton = ROOT.TLorentzVector(0, 0, 41, np.sqrt(41**2 + 0.938**2))  # Proton with 41 GeV in z-direction
initial_electron = ROOT.TLorentzVector(0, 0, -5, np.sqrt(5**2 + 0.000511**2))  # Electron with 5 GeV in -z direction

# Function to calculate pseudorapidity
def calculate_eta(p):
    theta = np.arccos(p.Pz() / p.P())
    return -np.log(np.tan(theta / 2))

# Function to create a TLorentzVector from a TParticle
def create_lorentz_vector(particle):
    return ROOT.TLorentzVector(particle.Px(), particle.Py(), particle.Pz(), particle.Energy())

# Function to calculate kinematic variables
def calculate_kinematics(electron):
    q = initial_electron - electron
    Q2 = -q.M2()
    W = (initial_proton + q).M()
    x = Q2 / (2 * initial_proton.Dot(q))
    return W, Q2, x

# Loop over the entries in the tree
for event in tree:
    particles = event.particles

    # Get the relevant particles
    electron = create_lorentz_vector(particles[4])
    phi = create_lorentz_vector(particles[5])
    proton = create_lorentz_vector(particles[6])
    kaon1 = create_lorentz_vector(particles[7])
    kaon2 = create_lorentz_vector(particles[8])

    # Calculate kinematic variables
    W, Q2, x = calculate_kinematics(electron)

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

    # Append momenta
    electron_momenta.append(electron.P())
    phi_momenta.append(phi.P())
    proton_momenta.append(proton.P())
    kaon1_momenta.append(kaon1.P())
    kaon2_momenta.append(kaon2.P())

    # Append pseudorapidity
    electron_eta.append(calculate_eta(electron))
    phi_eta.append(calculate_eta(phi))
    proton_eta.append(calculate_eta(proton))
    kaon1_eta.append(calculate_eta(kaon1))
    kaon2_eta.append(calculate_eta(kaon2))

# Define binning parameters
momentum_bins = np.linspace(0, 41, 50)  # Example range for momenta
eta_bins = np.linspace(-2.5, 5, 200)       # Example range for pseudorapidity

# Plotting
plt.figure(figsize=(12, 12))

plt.subplot(3, 2, 1)
plt.hist(electron_momenta, bins=momentum_bins, alpha=0.5, label='Electron')
plt.hist(phi_momenta, bins=momentum_bins, alpha=0.5, label='Phi Meson')
plt.hist(proton_momenta, bins=momentum_bins, alpha=0.5, label='Proton')
plt.hist(kaon1_momenta, bins=momentum_bins, alpha=0.5, label='Kaon1')
plt.hist(kaon2_momenta, bins=momentum_bins, alpha=0.5, label='Kaon2')
plt.xlabel('Momentum (GeV)')
plt.ylabel('Counts')
plt.legend()
plt.title('Momentum Distribution')

plt.subplot(3, 2, 2)
plt.hist(electron_eta, bins=eta_bins, alpha=0.5, label='Electron')
plt.hist(phi_eta, bins=eta_bins, alpha=0.5, label='Phi Meson')
plt.hist(proton_eta, bins=eta_bins, alpha=0.5, label='Proton')
plt.hist(kaon1_eta, bins=eta_bins, alpha=0.5, label='Kaon1')
plt.hist(kaon2_eta, bins=eta_bins, alpha=0.5, label='Kaon2')
plt.xlabel('Pseudorapidity (η)')
plt.ylabel('Counts')
plt.legend()
plt.title('Pseudorapidity Distribution')

plt.subplot(3, 2, 3)
plt.hist(W_values, bins=50, alpha=0.5)
plt.xlabel('W (GeV)')
plt.ylabel('Counts')
plt.title('W Distribution')

plt.subplot(3, 2, 4)
plt.hist(Q2_values, bins=50, alpha=0.5)
plt.xlabel('Q² (GeV²)')
plt.ylabel('Counts')
plt.title('Q² Distribution')

plt.subplot(3, 2, 5)
plt.hist(x_values, bins=50, alpha=0.5)
plt.xlabel('x')
plt.ylabel('Counts')
plt.title('x Distribution')

plt.tight_layout()
plt.show()
