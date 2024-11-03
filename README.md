# MAVD
Differential Abundance Analysis  is a vital statistical task in the analysis of microbiome data. Addressing the distinctive challenges inherent in microbiome abundance data, such as compositional effects, zero inflation, overdispersion, and high dimensionality, this study introduces the Microbial Abundance Data Vector Decomposition Model (MAVD)

The main analytical steps are as follows:
1)	Data standardization or normalization preprocessing.
2)	Linear modeling of normal group data using the flat method.
3)	Determination of the number of principal components using Wold invariance.
4)	Extraction of the first j principal components of normal data to obtain the normal plane space N.
5)	Leave-one-out cross-validation to test reliability and obtain the residual vector L1.Nmat of normal data.
6)	Projection of case group matrices to obtain D.Dmat and N.Dmat.
7)	Definition of disease threshold, such as mean(|D.Dmat|) > 0.95 quantile_probs(|L1.Nmat|).
8)	Filtering and removal of vectors below the threshold.
