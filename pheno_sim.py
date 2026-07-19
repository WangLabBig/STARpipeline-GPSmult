"""
To generate simulated phenotypes with realistic genetic architectures, we constructed individual-level BMI and T2D phenotypes based on empirical PRS. 
The observed BMI-PRS and T2D-PRS values were standardized to have a mean of zero and a variance of one before simulation. 
For BMI, we assumed a SNP-based heritability ((h2)) of 0.50 and specified that the empirical BMI-PRS explained 18% of the total phenotypic variance ((R2 of PRS=0.18)). 
The PRS contribution was modeled as: GPRS=PRSBMI×(R2 of PRS)1/2, so that the variance explained by the PRS component was equal to the predefined (R2 of PRS). 
The remaining genetic component not captured by the PRS was simulated independently from a normal distribution with variance: Var(Gresidual)=h2-(R2 of PRS), representing additional genetic effects from variants not included or captured by the PRS. 
Environmental and unexplained residual effects were simulated independently with variance: 1-h2. 
The latent BMI liability was then generated as the sum of the PRS-derived genetic component, residual genetic component, and environmental component: LBMI=GPRS+Gresidual+E.
The resulting continuous liability score was transformed to an observed BMI scale using population parameters from European adults (mean BMI = 27.3 kg/m2 and standard deviation = 4.75 kg/m2): BMI=LBMI×4.75+27.3.

For T2D, we applied a similar liability-based framework. We assumed a total heritability of 0.50 and a PRS-explained variance of 0.10 (R2 of PRS=0.10). 
The residual genetic component was simulated with variance: Var(Gresidual)=h2-(R2 of PRS), and non-genetic variance (1-h2), respectively. 
To convert the continuous liability score into a binary disease phenotype, we applied a liability-threshold model. 
The resulting liability score was converted into binary disease status using a liability-threshold model calibrated to a population prevalence of 12%. 
Individuals with liability scores exceeding this threshold were classified as T2D cases, while all remaining individuals were classified as controls. 
"""

import numpy as np
from scipy import stats
import pandas as pd

# compute PRS for BMI and T2D, and rename them as BMI-score and T2D-score
df = pd.read_csv('simu.score', sep='\t')
prs_bmi = np.array(df['BMI-score'].tolist())
prs_t2d = np.array(df['T2D-score'].tolist())

np.random.seed(42)
n = 50000

r2_bmi, r2_t2d = 0.18, 0.1
h2_bmi, h2_t2d = 0.5, 0.5
k_t2d = 0.12

prs_bmi_s = (prs_bmi - prs_bmi.mean()) / prs_bmi.std()
prs_t2d_s = (prs_t2d - prs_t2d.mean()) / prs_t2d.std()

var_g_resid_bmi = h2_bmi - r2_bmi   # 0.60 - 0.07 = 0.53
var_g_resid_t2d = h2_t2d - r2_t2d   # 0.50 - 0.10 = 0.40

g_resid_bmi = np.random.normal(0, np.sqrt(var_g_resid_bmi), n)
g_resid_t2d = np.random.normal(0, np.sqrt(var_g_resid_t2d), n)

var_e_bmi = 1 - h2_bmi   # 0.40
var_e_t2d = 1 - h2_t2d   # 0.50

e_bmi = np.random.normal(0, np.sqrt(var_e_bmi), n)
e_t2d = np.random.normal(0, np.sqrt(var_e_t2d), n)

liability_bmi = prs_bmi_s * np.sqrt(r2_bmi) + g_resid_bmi + e_bmi
liability_t2d = prs_t2d_s * np.sqrt(r2_t2d) + g_resid_t2d + e_t2d

bmi_mean, bmi_sd = 27.3, 4.75
bmi = liability_bmi * bmi_sd + bmi_mean

threshold = stats.norm.ppf(1 - k_t2d)
t2d = (liability_t2d > threshold).astype(int)
df['T2D'] = t2d
df['BMI'] = bmi

print("===BMI===")
print(f"Mean: {bmi.mean():.2f}, SD: {bmi.std():.2f}")
print(f"PRS R2 in simulated data: {np.corrcoef(prs_bmi_s, liability_bmi)[0,1]**2:.3f}")

print("\n===T2D===")
print(f"Case fraction: {t2d.mean():.3f}  (target: {k_t2d})")
corr_cases = np.corrcoef(prs_t2d_s[t2d==1], liability_t2d[t2d==1])[0,1]**2
print(f"PRS R2 in cases (liability): {corr_cases:.3f}")

df[['IID', 'BMI', 'T2D']].to_csv('simulate.txt', sep='\t', index=False)
