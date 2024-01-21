# GaiaXPcorrection

GaiaDR3XPspectracorrectionV1 is an astronomy Python package to provide a correction
for the Gaia DR3 XP spectra. The corrected spectra can significantly eliminate the 
systematic errors of the original Gaia DR3 XP spectra, please refer to our 
paper (Huang et al. 2024) for more details. 

We can give robust corrections for sources in the roughly following range:
−0.5 < BP − RP < 2
3 < G < 17.5 
E(B − V) < 0.8
Sources that do not fall into this range are not tested for the correctness of the 
corrections given by the correction package.

More information for this package in
https://github.com/HiromonGON/GaiaXPcorrection/blob/main/GaiaDR3XPspectracorrectionV1/readme

Chinese Version in 
https://github.com/HiromonGON/GaiaXPcorrection/blob/main/GaiaDR3XPspectracorrectionV1/%E8%87%AA%E8%BF%B0%E6%96%87%E4%BB%B6

# How to install
From source 

This package can be installed from the source code after downloading it from the git repo: 
https://github.com/HiromonGON/GaiaXPcorrection

or download the source code from the permalink:
https://doi.org/10.12149/101376

# Quick start
Correction for one spectrum

from GaiaDR3XPspectracorrectionV1 import Gaia_Correction_V1
[flux_out,caution,C2,C3] = Gaia_Correction_V1.correction(flux_origin,G,error,Truncation=False,absolute_correction=True)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Input              Description           Format                     Default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
flux_origin        original spectrum     nparray                     none
G                  G magnitude           float                       none
error              flux error            nparray                     False
Truncation         Optional in GaiaXPy   bool                        False
absolute_correction absolute_correction  bool                        True
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Output             Description                          Format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
flux_out           corrected spectrum                 float64 nparray          
caution            Reliability of corrections         int64 nparray
C2                 derived from original spectrum     float64
C3                 derived from original spectrum     float64
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

