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
