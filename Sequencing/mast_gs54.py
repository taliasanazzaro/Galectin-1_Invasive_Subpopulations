import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import binom, false_discovery_control
from rpy2.robjects import pandas2ri, r, globalenv, conversion
import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

r('library(MAST)')
r('library(base)')

ad = sc.read_h5ad('GS54_cleaned.h5ad')

ad.obs['substrate'] = ['soft' if c.startswith('Soft') else 'stiff' for c in ad.obs['sample']]
ad.obs['gel_side'] = ['bottom' if c.endswith('Mi') else 'top' for c in ad.obs['sample']]

def mast(adata, test_obs, ref, test, cov_obs):
    random_effects = '\n'.join([f'{obs} <- factor(colData(sca)${obs})\ncolData(sca)${obs} <- {obs}' for obs in cov_obs])
    formula = '~'+test_obs+''.join([' + (1 | '+obs+')' for obs in cov_obs])
    if len(cov_obs) == 0:
        method = 'glm'
    else:
        method = 'glmer'

    # R function
    r(f'''
        find_de_MAST_RE <- function(adata_) {{
            sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay")

            # store the columns that we are interested in as factors
            {test_obs} <- factor(colData(sca)${test_obs})

            # set the reference level
            {test_obs} <- relevel({test_obs},'{ref}')
            colData(sca)${test_obs} <- {test_obs}

            # same for donors (which we need to model random effects)
            # formatted random_effects string to look like:
            # obs <- factor(colData(sca)$obs)
            # colData(sca)$obs <- obs

            {random_effects}

            # define and fit the model
            zlmCond <- zlm(formula = {formula}, 
                        sca=sca, 
                        method='{method}', 
                        ebayes=F)#, 
                        # strictConvergence=F,
                        #cfitArgsD=list(nAGQ = 0)) # to speed up calculations
            
            # perform likelihood-ratio test for the condition that we are interested in    
            summaryCond <- summary(zlmCond, doLRT='{test_obs}{test}')

            # get the table with log-fold changes and p-values
            summaryDt <- summaryCond$datatable
            res <- merge(summaryDt[contrast=='{test_obs}{test}' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                            summaryDt[contrast=='{test_obs}{test}' & component=='logFC', .(primerid, coef)],
                            by='primerid') # logFC coefficients
            # MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
            res[,coef:=res[,coef]/log(2)]
            # do multiple testing correction
            res[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
            # res = res[res$FDR<0.05,, drop=F]

            return(res)
        }}
    ''')

    df = pd.DataFrame(adata.X.A, index=adata.obs_names, columns=adata.var_names)
    df = df.join(adata.obs[[test_obs]+cov_obs])
    ad_r = sc.AnnData(df[adata.var_names], obs=df.drop(columns=adata.var_names).astype('string'))

    globalenv['adata'] = ad_r

    return r('res <- find_de_MAST_RE(adata)').set_index('primerid')



min_cell_frac = 0.1

ad_sub = ad[:, (ad.X != 0).mean(axis=0) > 0.1]

i=1
n=6

print(f'{i}/{n}')
#mast(ad_sub, 'substrate', 'soft', 'stiff', ['gel_side']).to_csv('mast_gs54/stiff_vs_soft.csv')
i += 1
print(f'{i}/{n}')
#mast(ad_sub, 'gel_side', 'top', 'bottom', ['substrate']).to_csv('mast_gs54/top_vs_bottom.csv')

i += 1
print(f'{i}/{n}')
mast(ad_sub[ad_sub.obs['substrate'] == 'stiff'], 'sample', 'StiffNon', 'StiffMi', []).to_csv('mast_gs54/stiffmi_vs_stiffnon.csv')
i += 1
print(f'{i}/{n}')
mast(ad_sub[ad_sub.obs['substrate'] == 'soft'], 'sample', 'SoftNon', 'SoftMi', []).to_csv('mast_gs54/softmi_vs_softnon.csv')

i += 1
print(f'{i}/{n}')
mast(ad_sub[ad_sub.obs['gel_side'] == 'bottom'], 'sample', 'SoftMi', 'StiffMi', []).to_csv('mast_gs54/stiffmi_vs_softmi.csv')
i += 1
print(f'{i}/{n}')
mast(ad_sub[ad_sub.obs['gel_side'] == 'top'], 'sample', 'SoftNon', 'StiffNon', []).to_csv('mast_gs54/stiffnon_vs_softnon.csv')