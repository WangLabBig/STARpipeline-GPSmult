### Simulation of genotypes based on Hapnest
- Details access: https://github.com/intervene-EU-H2020/synthetic_data
- configure file is "config.yaml"


### Simulation of phenotypes based on polygenic risk scores

code: pheno_sim.py

To generate simulated phenotypes with realistic genetic architectures, we constructed individual-level BMI and type 2 diabetes (T2D) phenotypes based on empirical polygenic risk scores (PRS). The observed BMI-PRS and T2D-PRS values were standardized to have a mean of zero and a variance of one before simulation.

For BMI, we assumed a SNP-based heritability ($h^2$) of 0.50 and specified that the empirical BMI-PRS explained 18% of the total phenotypic variance ($R^2_{PRS}=0.18$). The PRS contribution was modeled as:


$$G_{PRS}=PRS_{BMI}\times\sqrt{R^2_{PRS}}$$


so that the variance explained by the PRS component was equal to the predefined ($R^2_{PRS}$). The remaining genetic component not captured by the PRS was simulated independently from a normal distribution with variance:

$$
Var(G_{residual})=h^2-R^2_{PRS}=0.50-0.18=0.32
$$

representing additional genetic effects from variants not included or captured by the PRS. Environmental and unexplained residual effects were simulated independently with variance:

$$
Var(E)=1-h^2=0.50
$$

The latent BMI liability was then generated as the sum of the PRS-derived genetic component, residual genetic component, and environmental component:

$$
L_{BMI}=G_{PRS}+G_{residual}+E
$$

The resulting continuous liability score was transformed to an observed BMI scale using population parameters from European adults (mean BMI = 27.3 kg/$m^2$ and standard deviation = 4.75 kg/$m^2$):

$$
BMI=L_{BMI}\times4.75+27.3
$$

For T2D, we applied a similar liability-based framework. We assumed a total heritability of 0.50 and a PRS-explained variance of 0.10 ($R^2_{PRS}=0.10$). The residual genetic component was simulated with variance:

$$
Var(G_{residual})=0.50-0.10=0.40
$$

and the environmental component was generated with variance:

$$
Var(E)=1-0.50=0.50
$$

The latent T2D liability was calculated as:

$$
L_{T2D}=PRS_{T2D}\times\sqrt{R^2_{PRS}}+G_{residual}+E
$$

To convert the continuous liability score into a binary disease phenotype, we applied a liability-threshold model. Assuming a population prevalence of 12% for T2D in individuals of European ancestry, the disease threshold was defined as the 88th percentile of the standard normal distribution:

$$
Threshold=\Phi^{-1}(1-0.12)
$$

Individuals with liability scores exceeding this threshold were classified as T2D cases, while all remaining individuals were classified as controls.

This simulation framework generated phenotypes in which the empirical PRS accounted for a predefined proportion of phenotypic variance, while additional genetic and environmental components were incorporated to achieve realistic heritability and disease prevalence distributions.
