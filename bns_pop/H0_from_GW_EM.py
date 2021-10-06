#!/usr/bin/env python
from gwcosmo.utilities.standard_cosmology import fast_cosmology
from astropy.table import Table
from gwpy.table import Table as gwTable
from ligo.skymap.distance import marginal_pdf
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import quad
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as constants
from scipy.stats import gaussian_kde
import pkg_resources
import pickle
from gwcosmo.utilities.standard_cosmology import redshift_prior
from tqdm import tqdm
from scipy.stats import ncx2, norm, truncnorm


def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[max(y), mean, sigma])
    return popt
def gauss(x, A, x0, sigma):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
    
def dL_by_z_H0(z,H0,Om0):
    speed_of_light = constants.c.to('km/s').value
    cosmo = fast_cosmology(Omega_m=Om0)
    return cosmo.dl_zH0(z, H0)/(1+z) + speed_of_light*(1+z)/(H0*E(z,Om0))

def E(z,Om):
    return np.sqrt(Om*(1+z)**3 + (1.0-Om))

# Define the H0 for the analysis
H0_array = np.linspace(20,200,100) 
H0_true = 70
Om0 = 0.308
sigma_z = 1e-4
# Redshift at which you want to evaluate your det. probability
z_array = np.linspace(0, 1,     100000)
d_array = np.linspace(0, 10000, 100000)
data_path = pkg_resources.resource_filename('gwcosmo', 'data/')
cosmo_fast = fast_cosmology(Omega_m=Om0)

# detection probability
pdet = pickle.load(open(data_path+'O2PSD_BNS_Nsamps20000_full_waveform_snr_12.0.p', 'rb'))
p_z = redshift_prior()
posterior_total = np.ones_like(H0_array)
posteriors = []
table = gwTable.read("coinc_10000.xml", format="ligolw", tablename="sim_inspiral")

for k in tqdm(range(442)):
    # print(k)
    # Distance posterior from BAYESTAR skymap
    m = Table.read("fits/"+str(k)+".fits", format="fits")
    prob = m["PROBDENSITY"]
    mu = m["DISTMU"]
    sigma = m["DISTSIGMA"]
    norm = m["DISTNORM"]

    pdf = marginal_pdf(d_array, prob, mu, sigma, norm)
    A, mean_d, sigma_d = gauss_fit(d_array, pdf)

    # p_dist_fit = gauss(d_array, A, mean_d, sigma_d)
    # p_dist_fit/=np.trapz(p_dist_fit,d_array)
    # p_dist_sample = np.random.normal(mean_d, sigma_d, 1000)
    a = (0.0 - mean_d) / sigma_d # boundary so samples don't go below 0
    p_dist_sample = truncnorm.rvs(a, 10000, loc=mean_d, scale=sigma_d, size=1000)
    
    # EM data
    cosmo = FlatLambdaCDM(H0=H0_true, Om0=Om0)
    d_true = table[k]['distance']
    # d_true = mean_d
    z_true = z_at_value(cosmo.luminosity_distance, d_true*u.Mpc)
    a = (0.0 - z_true) / sigma_z # boundary so samples don't go below 0
    z_galaxy = truncnorm.rvs(a, 5, loc=z_true, scale=sigma_z, size=1000)
    
    # Redshift at which you want to evaluate your det. probability
    posterior = np.ones_like(H0_array)
    for i,H0 in tqdm(enumerate(H0_array)):
        beta = np.trapz(pdet.pD_zH0_eval(z_array,H0)*p_z.p_z(z_array),z_array)
    
        # from distance to redshift distribution
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        dLs = cosmo.luminosity_distance(z_array).to(u.Mpc).value
        print(min(dLs), max(dLs))
        print(min(p_dist_sample), max(p_dist_sample))

        z_at_dL = interp1d(dLs, z_array)
        redshift = z_at_dL(p_dist_sample)
        # include dz/ddL jacobian and remove dl^2 prior from the posterior
        weights = 1/(dL_by_z_H0(redshift,H0,Om0)*cosmo_fast.dl_zH0(redshift,H0)**2)
        ww = np.sum(weights)
        z_likelihood = gaussian_kde(redshift,weights=weights)
        # marginalize over redshift
        posterior[i]=np.sum(z_likelihood(z_galaxy)*p_z.p_z(z_galaxy))*ww/beta

    # Normalize
    posterior/=np.trapz(posterior,H0_array)
    fig = plt.figure()
    plt.plot(H0_array,posterior, label=str(mean_d))
    plt.xlabel('r$H_0$ [km/Mpc/s]')
    plt.axvline(70)
    plt.ylabel('pdf')
    plt.legend()
    plt.savefig(str(k)+".png")
    
    posteriors.append(posterior)
    posterior_total *= posterior

fig = plt.figure()
for k in tqdm(range(442)):
    plt.plot(H0_array,posteriors[i])
plt.xlabel('r$H_0$ [km/Mpc/s]')
plt.axvline(70)
plt.ylabel('pdf')
plt.savefig("single.png")

fig = plt.figure()
plt.plot(H0_array,posterior_total)
plt.xlabel('r$H_0$ [km/Mpc/s]')
plt.axvline(70)
plt.ylabel('pdf')

plt.legend()
plt.savefig("H0_combined.png")