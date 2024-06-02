# Introduction
This project presents the analysis process and results of the SpatialSPRITE article (https://doi.org/10.1101/2024.05.07.592900).

Citation: Yuang Ma, Bo Gou, Yuetong Xu, Muya Shu, Falong Lu, Xiang Li. Embryo spatial 3D genomics. biorxiv.

Note:
1. Chinese version of the documentation: README-CH.md
2. The file paths in the scripts have not been modified; please modify them before running.

<br>

# Result
## 1. Distinguishing Sample Spots and Background Spots on the Slice
### 1.1 Determining Thresholds for Sample Spots and Background Spots
We use the end num within a spot as the thresholds to classify the spots on the slice into sample or background. Interaction information within background spots will be considered noise and filtered out, not included in subsequent analyses. To ensure accurate classification, we compared the classification results using different end num values as thresholds with the fluorescence results (Figure 1). Each point in Figure 1 represents a spot, where blue points represent spots with end num < threshold, classified as background, and red points represent spots with end num > threshold, classified as sample. Ultimately, the threshold for Sample 1 (ultrasound 10 minutes) is end num = 4000, and for Sample 2 (ultrasound 6 minutes) is end num = 5000.

```
script:
Distinguish_SampleBack/Distinguish_SampleBack.py
input:
Data/clusters_MS0612-5
Data/clusters_MS0612-3.odd70
output:
Data/clusters_MS0612-5.clean.sprite
Data/clusters_MS0612-3.odd70.clean.sprite
Distinguish_SampleBack/clusters_MS0612-5.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.pdf
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.pdf
Distinguish_SampleBack/clusters_MS0612-5.sample.spotid
Distinguish_SampleBack/clusters_MS0612-3.odd70.sample.spotid
Data/clusters_MS0612-5.sample.sprite
Data/clusters_MS0612-3.odd70.sample.sprite
```

<img src="Doc/Fig1_SampleBack.png" alt="fig 1" />
<center>Fig. 1. Classification results of spots under different thresholds</center><br>


### 1.2 Statistical Analysis of Spot Data
After filtering, some non-sample region spots still had an end num > threshold. We manually deleted these spots by comparing them with the fluorescence results. For details, refer to Distinguish_SampleBack.py. Ultimately, there are 3106 spots classified as sample and 1694 spots classified as background in the Sample 1 slice; in the Sample 2 slice, 1853 spots are classified as sample, 1647 spots are classified as background. Additionally, we also counted the cluster num and end num within each spot for both samples (Figure 2). In Figure 2, the spots are sorted according to either cluster num or end num.

```
script:
Distinguish_SampleBack/Heatmap_ClusterNum_EndNum.py
Distinguish_SampleBack/Scatter_ClusterNum_EndNum.py
input:
Distinguish_SampleBack/clusters_MS0612-5.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.tsv
output:
Heatmap_ClusterNum_EndNum.10min.pdf
Heatmap_ClusterNum_EndNum.6min.pdf
Scatter_ClusterNum_EndNum.pdf
```

<img src="Doc/Fig2_ClusterNum_EndNum_Spot.png" alt="fig 2" />
<center>Fig. 2. Cluster num and end num for each spot in the slice</center>
<p><br>


### 1.3 Proportion of Valid Information in Sample Region Spots
We consider the interaction information within background region spots as noise and calculate the proportion of valid information in sample region spots by '(sample-back)/sample'. Specifically, we first calculate the average end num within background region spots to represent background noise. Then, we calculate the proportion of valid end num within sample region spots. We display the proportion of valid information in all sample region spots in the slice using a histogram (Figure 3).

```
script:
Distinguish_SampleBack/Histogram_ValidInformRatio.py
input:
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-5.sample.spotid
Distinguish_SampleBack/clusters_MS0612-3.odd70.sample.spotid
output:
Distinguish_SampleBack/Histogram_ValidInformRatio.pdf
```

<img src="Doc/Fig3_Histogram_ValidInformRatio.png" alt="fig 3" />
<center>Fig. 3. Distribution of Valid Information Proportion in Sample Region Spots</center><br>


## 2. SPRITE Data Filtering
Since this experiment does not consider interactions between chromosomes, we first classify the multi-chrom SPRITE (inter-chromosome) within the sample region spots into single-chrom SPRITE (intra-chromosome). Then, we filter out adjacent (30bp) ends within the clusters and large clusters. Filtering adjacent ends is to avoid experimental errors caused by random primers, and filtering large clusters is to eliminate chromosome fragments that are not sufficiently broken. Specifically, we divide the chromosomes into 1MB bins and then count the number of bins that contain ends from each cluster. If this number exceeds 40% of the total number of bins, the cluster is considered to represent a large chromatin fragment that is not sufficiently broken and will be filtered out.

```
script:
Data_Filtering/Filter_ClosedEnds_BigCluster.py
input:
Data/clusters_MS0612-5.sample.sprite
Data/clusters_MS0612-3.odd70.sample.sprite
output:
Data/clusters_MS0612-5.sample.intra.sprite
Data/clusters_MS0612-3.odd70.sample.intra.sprite
Data/clusters_MS0612-5.sample.intra.FRI.F1.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.sprite
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
```

<center>Table. 1. The Data of SPRITE at Different Stages (Cluster Num)</center>

|       | All   | Sample    | Intra | F1    | FB (40%)  |
| ----  | ----  | ----      | ----  | ----  | ----:     |
|Sample 1| 2,562,880 | 2,248,559 | 1,829,574 | 1,335,255 | 1,138,141 |
|Sample 2| 2,266,648 | 1,806,029 | 1,486,294 | 1,049,433 | 932,713 |

<p><br>


## 3. Drawing Interaction Heatmaps of Chromosomes
Firstly, we convert the filtered SPRITE data from multi-ends interaction format to double-ends format using a combination method. The converted data volume is shown in Table 2. Then, we utilize Juicer Tools Pre and hic2cool to convert the interaction data into cool format (interaction matrix). Subsequently, we normalize the matrix using cooler and draw interaction heatmaps for each chromosome (Figure 4).

```
script:
Contact_Heatmap/Prepare_Data.py
Contact_Heatmap/Heatmap_Contact.py
input:
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
output:
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.pairs
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.hic
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/Heatmap_Contact_MS0612-5.pdf
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.pairs
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.hic
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/Heatmap_Contact_MS0612-3.odd70.pdf
```

<center>Table. 2. Interaction number of Samples in Different Formats</center>

|      | cluster num（multi-ends）| contact num（double-ends）|
| ---- | ----:  | ----: | 
|Sample 1| 1,138,141 | 815,239,083 |
|Sample 2|   932,713 | 430,938,706 |

<p><br>

<img src="Doc/Fig4_Contact_Heatmap.png" alt="fig 4" />
<center>Fig. 4. Interaction Heatmaps of Chromosomes 4 and 6</center><br>


## 4. Calculating A/B Compartment
We use the data in cool format as input and perform PCA on the chromosomes at a 1Mbp resolution using the eigdecomp.cis_eig function from cooltools. The eigenvectors are sorted based on their similarity (Pearson correlation) with GC content. The eigenvector with the highest similarity is identified as E1. Subsequently, chromosome regions with E1 values greater than 0 are classified as A compartments, and those with E1 values less than 0 are classified as B compartments (Figure 5).

```
script:
AB_Compartment/Plot_ABCompartment.py
input:
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
output:
AB_Compartment/Lineplot_AB_Compartment_MS0612-5.pdf
AB_Compartment/Lineplot_AB_Compartment_MS0612-3.odd70.pdf
```

<img src="Doc/Fig5_AB_Compartment.png" alt="fig 5" />
<center>Fig. 5. A/B Compartments of Chromosomes 4 and 6</center><br>


## 5. Clustering Spots
We perform k-means clustering on the spots in the sample regions using the number of interactions (end num) within chromosomal intervals as features.

### 5.1 Constructing the Spot-Bin Interaction Matrix
We constructed a matrix (Sample 1: 3106 * 2737, Sample 2: 1853 * 2737) with rows representing spot IDs and columns representing bin IDs, counting the end num for each chromosomal interval (bin, 1 MB) within the spots. Then, we calculated the missing rates for spots and bins (Figure 6). We filtered out spots and bins with high missing rates (Top 5%) (Figure 7, Sample 1: 2950 * 2451, Sample 2: 1760 * 2450).

```
script:
KMeans/Filter_SpotBin_MissRate.py
input:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
output:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.matrix
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/Histplot_MS0612-5_FB_0.4.missrate.pdf
KMeans/Histplot_MS0612-5_FB_0.4.missrate.filter5%.pdf
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/Histplot_MS0612-3.odd70_FB_0.4.missrate.pdf
KMeans/Histplot_MS0612-3.odd70_FB_0.4.missrate.filter5%.pdf
```

<img src="Doc/Fig6_Histplot_Missrate.png" alt="fig 6" />
<center>Fig. 6. Missing Rates of Spots and Bins in the Spot-Bin Interaction Matrix</center><br>

<img src="Doc/Fig7_Histplot_FilterMissrate.png" alt="fig 7" />
<center>Fig. 7. Missing Rates of Spots and Bins in the Filtered Matrix (Top 5%)</center><br>


### 5.2 Normalization and Imputation
First, we clip the bins with extreme interaction counts (Top 0.5%) within spots. Then, we use z-score normalization, treating bins with no interactions (end num=0) as missing (None). Next, we impute data using interaction information from neighboring (N=2) spots (Figure 8). Finally, since the imputation process might alter the distribution of interaction data within spots, the script uses qnorm to re-normalize the spot-bin interaction matrix.

```
script:
KMeans/Normalized_Imputation.py
input:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.matrix
output:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
```

<center><img src="Doc/Fig8_Imputation.png" width = "200" height = "240" alt="fig 8" /></center>
<center>Fig. 8. Imputation of center spot using neighboring spots (N=2)</center><br>


### 5.3 PCA Dimensionality Reduction
We use PCA to reduce the dimensionality of the spot-bin matrix, decreasing the number of columns while increasing the amount of information per feature. We plotted the cumulative explained variance for each dimension after PCA (Figure 9). The results show a noticeable elbow at dimension=200, indicating that dimensions beyond 200 explain less variance. Therefore, we select the first 200 dimensions after PCA as the features for spots in subsequent clustering analysis, resulting in Sample 1: 2893 * 200, Sample 2: 1705 * 200. Since PCA does not support None values, we excluded spots that still had missing values after imputation, hence the number of spots is less than 3.1.

```
script:
KMeans/PCA.py
input:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
output:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.pdf
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.pdf
```

<center><img src="Doc/Fig9_PCA.png" alt="fig 9" /></center>
<center>Fig. 9. Relationship between PCA dimensions and explained variance</center><br>


### 5.4 K-means Clustering
We input the dimension-reduced spot-bin matrix into the K-means method (K=30) for clustering spots. The clustering results were unstable with different seeds, so we selected the results whose spatial positions and shapes of clusters best matched the organ data (MOSTA data). After selection, the parameters for Sample 1 (MS0612-5) were K=30, seed=31; for Sample 2 (MS0612-3.odd70), K=30, seed=11. Since some clusters contained spatially non-adjacent spots, we split them into multiple classes. Thus, the final number of classes is greater than 30. Sample 1 is divided into 35 classes (Figure 10), and Sample 2 is divided into 32 classes (Figure 11). Additionally, we associated the clusters with organs (Figure 12). For Sample 2, no suitable reference annotation was found, so the figure legend contains the authors' speculation.

```
script:
KMeans/KMeans.py
KMeans/Scatter_KMeans_Result.py
input:
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
output:
KMeans/KMeans_MS0612-5.label
KMeans/KMeans_MS0612-5.png
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-5.split.png
KMeans/Scatter_KMeans_Result.MS0612-5.split.pdf
KMeans/KMeans_MS0612-3.odd70.label
KMeans/KMeans_MS0612-3.odd70.png
KMeans/KMeans_MS0612-3.odd70.split.label
KMeans/KMeans_MS0612-3.odd70.split.png
KMeans/Scatter_KMeans_Result.MS0612-3.odd70.split.pdf
```

<center><img src="Doc/Fig10_KMeans_MS0612-5.split.png" width = "400" height = "600" alt="fig 10" /></center>
<center>Fig. 10. K-means Clustering Results and Spatial Distribution of Sample 1</center><br>

<center><img src="Doc/Fig11_KMeans_MS0612-3.odd70.split.png" width = "400" height = "600" alt="fig 11" /></center>
<center>Fig. 11. K-means Clustering Results and Spatial Distribution of Sample 2</center><br>

<center><img src="Doc/Fig12_KMeans.png" alt="fig 12" /></center>
<center>Fig. 12. K-means Clustering Results and Corresponding Organs</center><br>


## 6. Computing Interaction Frequencies at Different Distances for Various Tissues
First, we extracted SPRITE data for each cluster category based on the K-means results and converted it into a cool file format. Then, utilizing the expected_cis function from cooltools, we computed the relationship between genomic-level interaction distance and interaction frequency for each tissue, with the liver showing particular distinctiveness (Fig. 13). The results revealed that different tissues exhibit distinct patterns of interaction frequency at various distances. We categorized interactions into short-range (≤ 10 MB), mid-range (10 MB ＜ x ≤ 40 MB), and long-range (＞ 40 MB). In short-range interactions, the abdominal region exhibited higher interaction frequencies, while facial interactions were moderate, and brain interactions were lower. In mid-range interactions, the liver showed higher interaction frequencies compared to other tissues. In long-range interactions, the brain and facial regions exhibited higher interaction frequencies, while the abdominal region showed moderate frequencies, and the liver exhibited lower frequencies.

```
script:
Contact_Distance/Prepare_TissueData.py
Contact_Distance/Prepare_CoolData.py
Contact_Distance/Lineplot_ContactDistance.py
Contact_Distance/Heatmap_ContactDistance.py
input:
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
output:
Data/Tissue_Pairs/*
Contact_Distance/Data_HiC/*
Contact_Distance/Data_Cool/*
Contact_Distance/Contact_Distance.MS0612-5.tsv
Contact_Distance/Contact_Distance.MS0612-3.odd70.tsv
Contact_Distance/Lineplot_ContactDistance.MS0612-5.pdf
Contact_Distance/Lineplot_ContactDistance.MS0612-3.odd70.pdf
Contact_Distance/Heatmap_ContactDistance.MS0612-5.pdf
Contact_Distance/Heatmap_ContactDistance.MS0612-3.odd70.pdf
```

<center><img src="Doc/Fig13_Contact_Distance.png" alt="fig 13" /></center>
<center>Fig. 13. The interaction frequency at different distances in various tissues</center><br>


## 7. Analyzing Compartmentalization Strength of Various Tissues
### 7.1 Calculating A/B Compartments for Various Tissues
Firstly, we extracted SPRITE data for each category based on the K-means results. Then, following previous steps, we calculated the E1 values and A/B compartments for each chromosome within each category. Subsequently, we computed the Pearson correlation coefficient between the E1 values of chromosomes within categories and those of embryos (Fig. 14). The results indicated significant differences in E1 values for some chromosomes between categories and embryos. Hence, we focused on the top 10 chromosomes (2, 4, 5, 6, 8, 9, 11, 15, 17, 19) with high correlation for further analysis.

```
script:
Compartmentalization_Strength/Prepare_CoolData.py
Compartmentalization_Strength/Plot_ABCompartment.py
input:
Data/Tissue_Pairs/*
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
output:
Compartmentalization_Strength/Data_HiC/*
Compartmentalization_Strength/Data_Cool/*
Compartmentalization_Strength/AB_Compartment/*
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.value
Compartmentalization_Strength/E1_Value/MS0612-5.e1.chrom.pearson
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.chrom.pearson
Compartmentalization_Strength/E1_Value/Boxplot_MS0612-5.e1.chrom.pearson.pdf
Compartmentalization_Strength/E1_Value/Boxplot_MS0612-3.odd70.e1.chrom.pearson.pdf
```

<center><img src="Doc/Fig14_E1_Pearson.png" alt="fig 14" /></center>
<center>Fig. 14. The correlation coefficient between the E1 values of each tissue and the embryo</center><br>


### 7.2 Calculating Compartmentalization Strength of Various Tissues
We employed the saddle function from cooltools to partition chromosomes of each tissue into 40 bins based on embryonic E1 values and calculated the ratio of interaction frequencies between bins and their expected values. The results were presented in saddle plots (Fig. 14A). Bins were sorted based on their E1 values along the x and y axes, where the four quadrants represented the interaction frequency ratios between BB, BA, AB, and AA compartments, respectively. Subsequently, we utilized the saddle_strength function from cooltools to compute the compartmentalization strength of each tissue (formula, Fig. 14B). Here, AA, BB, AB, and BA denoted the average interaction frequency ratios within the 5*5 regions corresponding to the four quadrants of the saddle plot (extend=5, Fig. 14A). The results indicated higher compartmentalization strength in the embryonic abdominal region, moderate strength in the facial region, and lower strength in the brain region.

$$Strength=\frac{𝐵𝐵+AA}{𝐵A+A𝐵}$$

```
script:
Compartmentalization_Strength/Compartment_Strength.py
input:
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.value
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
Compartmentalization_Strength/Data_Cool/*
output:
Compartmentalization_Strength/Saddle.Strength_MS0612-5.pdf
Compartmentalization_Strength/Saddle.Strength_MS0612-3.odd70.pdf
Compartmentalization_Strength/Strength_MS0612-5.tsv
Compartmentalization_Strength/Strength_MS0612-3.odd70.tsv
Compartmentalization_Strength/Scatter.Strength_MS0612-5.pdf
Compartmentalization_Strength/Scatter.Strength_MS0612-3.odd70.pdf
```

<center><img src="Doc/Fig15_Compartmentalization_Strength.png" alt="fig 15" /></center>
<center>Fig. 15. Saddle-plots and compartmentalization strengths of various tissues</center><br>


## 8. Analyzing Interaction Frequency Features of Tissues with Different Compartmentalization Strength
We selected 2, 3, and 5 tissues from Sample 1 for each category of high, medium, and low compartmentalization strength, respectively (Fig. 16A). Then, on a compartment-by-compartment basis, we calculated the interaction frequencies within and between compartments of tissues with different strengths, and performed bootstrapping tests. The results showed that with increasing compartmentalization strength, inter-compartmental interaction intensity decreased (Fig. 16BC), while intra-compartmental interaction intensity increased (Fig. 16DE).

```
script:
Compartment_Contact/Calculate_Compartment_Border.py
Compartment_Contact/Calculate_Compartment_ContactFrequency.py
Compartment_Contact/Violinplot_Compartment_ContactFrequency.py
input:
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.value
Data/Tissue_Pairs/*
output:
Compartment_Contact/Compartment_AB_border.MS0612-5.tsv
Compartment_Contact/Compartment_AB_border.MS0612-3.odd70.tsv
Data/Tissue_Chrom_Pairs/*
Compartment_Contact/Compartment_ContactNum_Intra.MS0612-5.tsv
Compartment_Contact/Compartment_ContactNum_Inter.MS0612-5.tsv
Compartment_Contact/Violinplot_CompartmentA_ContactFrequency_Inter.MS0612-5.pdf
Compartment_Contact/Violinplot_CompartmentB_ContactFrequency_Inter.MS0612-5.pdf
Compartment_Contact/Violinplot_CompartmentA_ContactFrequency_Intra.MS0612-5.pdf
Compartment_Contact/Violinplot_CompartmentB_ContactFrequency_Intra.MS0612-5.pdf
Compartment_Contact/Bootstrapping_CompartmentA_ContactFrequency_Inter.MS0612-5.pdf
Compartment_Contact/Bootstrapping_CompartmentB_ContactFrequency_Inter.MS0612-5.pdf
Compartment_Contact/Bootstrapping_CompartmentA_ContactFrequency_Intra.MS0612-5.pdf
Compartment_Contact/Bootstrapping_CompartmentB_ContactFrequency_Intra.MS0612-5.pdf
```

<center><img src="Doc/Fig16_Compartment_ContactFrequency.png" alt="fig 16" /></center>
<center>Fig. 16. The distribution of interaction frequencies among tissues with different compartmentalization strengths</center><br>


## 9. Analyzing Characteristics of Long-Distance, Inter-Compartment Interactions in Liver Tissue
Due to the significant decrease in long-distance interactions (> 40 MB) in liver tissue (Label ID = 1) compared to other tissues (Fig. 13A), we specifically analyzed the interaction differences between the 11 tissues of Sample 1 within the same type of compartment (Fig. 17). The colors in the graph correspond to interaction frequency, with darker lines indicating higher frequency. The results showed a significant reduction in B-B compartment interaction frequency in liver tissue on chromosomes 2, 5, and 9.

```
script:
Long_Compartment_Contact/Prepare_LongContact.py
Long_Compartment_Contact/Lineplot_LongContact.py
input:
Data/Tissue_Pairs/*
Data/Tissue_Chrom_Pairs/*
Compartment_Contact/Compartment_AB_border.MS0612-5.tsv
Long_Compartment_Contact/mm10_chrom_sizes.txt
output:
Data/Tissue_Chrom_Pairs_LongDistance/*
Long_Compartment_Contact/Compartment_ContactNum_Inter.LongDistance.MS0612-5.tsv
Long_Compartment_Contact/Compartment_ContactNum_Inter.LongDistance.MS0612-5.pdf
```

<center><img src="Doc/Fig17_LongDistance_Inter-Compartment_Contact.png" alt="fig 17" /></center>
<center>Fig. 17. The differences in long-distance, inter-compartment interaction frequencies among different tissues</center><br>


## 10. Computing the 3D Structure of Chromosomes in Various Tissues
We used the Hickit software to calculate the three-dimensional structure of 10 chromosomes across 11 tissues at a resolution of 1 MB (Fig. 18). Additionally, we labeled the compartment categories (A/B) to which chromosome segments belonged.

```
script:
3D_Genome/Prepare_Data.py
3D_Genome/Scatter_3DGenome.py
3D_Genome/Scatter_3DGenome_AB.py
input:
Data/Tissue_Chrom_Pairs/*
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
output:
Data/Tissue_Chrom_Pairs_3DGenome/*
3D_Genome/Data_Hickit/*
3D_Genome/Scatter_3D_Structure/*
3D_Genome/Scatter_3D_Structure_AB/*
```

<center><img src="Doc/Fig18_Liver_Chrom8_3DStructure.png" alt="fig 18" /></center>
<center>Fig. 18. The 3D structure of chromosome 8 in liver tissue (Label ID=1)</center><br>