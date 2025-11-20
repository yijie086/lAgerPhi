import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file
file = ROOT.TFile("outputs/solid-ep-phi.ep-phi.composite.run01001-lumi100.root")

# Get the TTree
tree = file.Get("lAger")

# Initialize lists to store momenta and angles
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

# Function to calculate theta in degrees
def calculate_theta(p):
    return np.degrees(np.arccos(p.Pz() / p.P()))

# Loop over the entries in the tree
for event in tree:
    particles = event.particles
    
    # Get the relevant particles
    electron = particles[4]
    phi = particles[5]
    proton = particles[6]
    kaon1 = particles[7]
    kaon2 = particles[8]
    
    # Append momenta
    electron_momenta.append(electron.P())
    phi_momenta.append(phi.P())
    proton_momenta.append(proton.P())
    kaon1_momenta.append(kaon1.P())
    kaon2_momenta.append(kaon2.P())
    
    # Append theta
    electron_theta.append(calculate_theta(electron))
    phi_theta.append(calculate_theta(phi))
    proton_theta.append(calculate_theta(proton))
    kaon1_theta.append(calculate_theta(kaon1))
    kaon2_theta.append(calculate_theta(kaon2))

# Define binning parameters
momentum_bins = np.linspace(0, 15, 50)  # Example range for momenta
theta_bins = np.linspace(0, 180, 50)    # Example range for theta in degrees

# Plotting
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.hist(electron_momenta, bins=momentum_bins, alpha=0.5, label='Electron')
plt.hist(phi_momenta, bins=momentum_bins, alpha=0.5, label='Phi Meson')
plt.hist(proton_momenta, bins=momentum_bins, alpha=0.5, label='Proton')
plt.hist(kaon1_momenta, bins=momentum_bins, alpha=0.5, label='Kaon1')
plt.hist(kaon2_momenta, bins=momentum_bins, alpha=0.5, label='Kaon2')
plt.xlabel('Momentum (GeV)')
plt.ylabel('Counts')
plt.legend()
plt.title('Momentum Distribution')

plt.subplot(1, 2, 2)
plt.hist(electron_theta, bins=theta_bins, alpha=0.5, label='Electron')
plt.hist(phi_theta, bins=theta_bins, alpha=0.5, label='Phi Meson')
plt.hist(proton_theta, bins=theta_bins, alpha=0.5, label='Proton')
plt.hist(kaon1_theta, bins=theta_bins, alpha=0.5, label='Kaon1')
plt.hist(kaon2_theta, bins=theta_bins, alpha=0.5, label='Kaon2')
plt.xlabel('Theta (degrees)')
plt.ylabel('Counts')
plt.legend()
plt.title('Theta Distribution')

plt.tight_layout()
plt.show()
