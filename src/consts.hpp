#ifndef consts_hpp
#define consts_hpp


#include "mathdefs.hpp"

namespace flexiblesusy {

// units

const double TeV = 1;
const double GeV = 1e-3 * TeV;
const double MeV = 1e-3 * GeV;

const double c_m_s = 299792458;	// speed of light in m/s
const double hbarc_MeVfm = 197.3269602;	// hbar c in MeV fm
const double ips = hbarc_MeVfm / (1e6*c_m_s) * GeV;
const double cm = 1e13 / (hbarc_MeVfm*MeV);
const double sec = c_m_s * 1e15 / (hbarc_MeVfm*MeV);
const double year = 365*24*60*60*sec;
const double Gyear = 1e9*year;

const double deg = pi / 180;

// from PDB 2000.

const double sinSqThetaW = 0.23117; // (at mZ, bar MS) 0.23117(16)
const double al = 1 / 137.035;
const double almZ = 1.0/129;
const double gYmZ = sqrt(4*pi * almZ/(1-sinSqThetaW));
const double g1mZ = sqrt(5.0/3) * gYmZ;
const double g2mZ = sqrt(4*pi * almZ/sinSqThetaW);
const double mW = 80.419  * GeV; // 80.419(56)  GeV
const double mZ = 91.1876 * GeV; // +- 0.0021 GeV
const double vv = sqrt(2*sqr(mZ)/(sqr(g2mZ)+sqr(gYmZ)));
const double GF = 1.16639e-5 / sqr(GeV); // 1.16639(1)e-5 GeV^-2
const double alsmt = 0.108979;		 // alphas5(167 GeV) from
				         // ../../b3s/pr/LR/mathe-common.h
const double g3mt = sqrt(4*pi * alsmt);
const double mt_MSbar = 167*GeV; // m_t^{bar MS} (m_t^{bar MS})

// Running quark masses at mW from Table 1 of hep-ph/0202055.
const double muMW      =   1.7	*MeV;
const double sigmamuMW =   0.5	*MeV;
const double mdMW      =   3.7	*MeV;
const double sigmamdMW =   0.9	*MeV;
const double msMW      =   0.072*GeV;
const double sigmamsMW =   0.017*GeV;
const double mcMW      =   0.667*GeV;
const double sigmamcMW =   0.027*GeV;
const double mbMW      =   2.906*GeV;
const double sigmambMW =   0.047*GeV;
const double mtMW      = 175.3  *GeV;
const double sigmamtMW =   5.1  *GeV;

// Lepton masses.  PDG2004.
const double me = 0.51099892*MeV;
const double sigmame = 0.00000004*MeV;
const double mmu = 105.658369*MeV;
const double sigmammu = 0.000009*MeV;
const double mtau = 1776.99*MeV;
const double sigmamtau = 0.275*MeV;

// MS bar quark masses at mW scale
// from Caravaglios, Roudeau, Stocchi, NPB633, 193 (2002), hep-ph/0202055
const double mb_MSbar_mW = 2.906*GeV; // +- 0.035 +- 0.031
const double ms_MSbar_mW = 0.072*GeV; // +- 0.017 +- 0.002
const double md_MSbar_mW = 3.7*MeV;   // +- 0.9

const double mBd = 5.279*GeV;
const double mBs = 5.370*GeV;

// from Becirevic et al., NPB634, 105 (2002), hep-ph/0112303
const double fBd_cent = 200*MeV;
const double fBd_err  = 30*MeV;
const double B1_cent = 0.87;
const double B2_cent = 0.82;
const double B3_cent = 1.02;
const double B4_cent = 1.16;
const double B5_cent = 1.91;

// Okamoto, hep-lat/0510113
const double BsBdscale = 1.21;	// +.047 -.035
const double fBs_cent = fBd_cent*BsBdscale;
const double B1s_cent = B1_cent;
const double B2s_cent = B2_cent;
const double B3s_cent = B3_cent;
const double B4s_cent = B4_cent;
const double B5s_cent = B5_cent;

//
const double B1K = 0.6;
const double B2K = 0.66;
const double B3K = 1.05;
const double B4K = 1.03;
const double B5K = 0.73;

const double FK  = 160*MeV; // Buchalla, Buras, Lautenbacher, RMP68, 1125 (96)
const double mK  = 497.648*MeV;	// +- 0.022*MeV, PDG2006

// Numbers from hep-ph/0112303
const double md = 5.4*MeV;	// RI-MOM mass at m_b
// const double ms = 130*MeV;
const double mc = 1.25*GeV;	// PDG 2000
// This value of ms may not be based on a reference.
const double ms = 100*MeV;	// RI-MOM mass at m_b
const double mb = 4.6*GeV;	// RI-MOM mass at m_b
const double mb_pole = 4.6*GeV;	// pole mass
const double mb_MSbar = 4.23*GeV; // m_b^{bar MS} (m_b^{bar MS})

const double ms_MSbar = 115*MeV; // +- 20 MeV, m_s(m_c), hep-ph/0306217
// This md was guessed from (2.9) and (3.2) in hep-ph/0306217.
const double md_MSbar = 6*MeV;	// m_d^{bar MS} (m_c)

} // namespace flexiblesusy

#endif // consts_hpp
