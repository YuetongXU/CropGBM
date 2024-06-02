## 介绍
本项目介绍了 SpatialSPRITE 文章（ https://doi.org/10.1101/2024.05.07.592900 ）的分析流程及结果。

**Citation**: Yuang Ma, Bo Gou, Yuetong Xu, Muya Shu, Falong Lu, Xiang Li. Embryo spatial 3D genomics. biorxiv.

**Note**: 
1. 英文版说明文档: README.md
2. 脚本中文件路径未修改，运行时请修改。

<br>


# 结果
## 1. 区分切片上的样本 spot 及背景 spot
### 1.1 确定样本 spot 及背景 spot 的阈值
我们以 spot 内 end num 为标准，划分切片上 spot 的类型（样本区 or 背景区）。其中，背景区 spot 内的互作信息将被视为噪音而过滤掉，不纳入后续的分析中。为了确保划分地准确，我们比较了不同 end num 作为阈值时的划分结果（图 1）与荧光结果的吻合程度。图 1 中每个点代表一个spot，其中蓝色点代表 end num &lt; threshold 的 spot，其被判定为背景区；红色点代表 end num &gt; threshold 的spot，其被判定为样本区。最终，样本1（超声 10min）的阈值为 end num=4000，样本2（超声 6min）的阈值为 end num=5000。

```
脚本：
Distinguish_SampleBack/Distinguish_SampleBack.py
输入：
Data/clusters_MS0612-5
Data/clusters_MS0612-3.odd70
输出：
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
<center>Fig. 1. 不同阈值下 spot 的划分结果</center><br>


### 1.2 统计 spot 的数据量
过滤后仍有部分非样本区 spot 的 end num &gt; threshold，这里我们对照荧光结果进行了手动删除，详情参见 Distinguish_SampleBack.py。最终，Sample 1 切片内有 3106 个 spot 隶属于样本区，1694 个 spot 隶属于背景区；Sample 2 切片内有 1853 个 spot 隶属于样本区，1647 个 spot 隶属于背景区。此外，我们还统计了两个样本各 spot 内的 cluster num 及 end num（图 2）。图 2 中我们将 spot 按照 cluster num 或 end num 排序。

```
脚本：
Distinguish_SampleBack/Heatmap_ClusterNum_EndNum.py
Distinguish_SampleBack/Scatter_ClusterNum_EndNum.py
输入：
Distinguish_SampleBack/clusters_MS0612-5.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.clusternum.tsv
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.tsv
输出：
Heatmap_ClusterNum_EndNum.10min.pdf
Heatmap_ClusterNum_EndNum.6min.pdf
Scatter_ClusterNum_EndNum.pdf
```

<img src="Doc/Fig2_ClusterNum_EndNum_Spot.png" alt="fig 2" />
<center>Fig. 2. 切片内各 spot 的 cluster num 及 end num 情况</center>
<p><br>


### 1.3 统计样本区 spot 内有效信息的比例
我们以背景区 spot 内的互作信息为噪音，统计了样本区 spot 内有效信息的比例 ((sample-back)/sample)。具体而言，我们首先计算背景区 spot 内 end num 的平均值，以此代表背景噪音，再计算样本区 spot 内有效 end num 的比例，然后我们将切片内所有样本区 spot 的有效信息比例以直方图的方式展示（图 3）。

```
脚本：
Distinguish_SampleBack/Histogram_ValidInformRatio.py
输入：
Distinguish_SampleBack/clusters_MS0612-5.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-3.odd70.spot.endnum.tsv
Distinguish_SampleBack/clusters_MS0612-5.sample.spotid
Distinguish_SampleBack/clusters_MS0612-3.odd70.sample.spotid
输出：
Distinguish_SampleBack/Histogram_ValidInformRatio.pdf
```

<img src="Doc/Fig3_Histogram_ValidInformRatio.png" alt="fig 3" />
<center>Fig. 3. 切片内样本区 spot 的有效信息比例分布</center><br>


## 2. SPRITE 数据过滤
因为本实验不考虑染色体之间的互作，所以我们首先将样本区 spot 内的 mutil-chrom SPRITE (inter-chromosome) 划分为 single-chrom SPRITE (intra-chromosome)。然后过滤 cluster 中的相邻 (30bp) end 和大 cluster。其中，过滤相邻 end 是为了避免随机引物引起的实验误差；过滤大 cluster 是为了剔除未充分断裂的染色体片段。具体来说，我们将染色体分成长度为 1MB 的 bin，然后统计 cluster 的 end 所在 bin 的数量。如果数量超过染色质总 bin 数的 40%，则认为该 cluster 代表一个未充分断裂的染色质大片段，将被过滤掉。

```
脚本：
Data_Filtering/Filter_ClosedEnds_BigCluster.py
输入：
Data/clusters_MS0612-5.sample.sprite
Data/clusters_MS0612-3.odd70.sample.sprite
输出：
Data/clusters_MS0612-5.sample.intra.sprite
Data/clusters_MS0612-3.odd70.sample.intra.sprite
Data/clusters_MS0612-5.sample.intra.FRI.F1.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.sprite
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
```

<center>Table. 1. 不同阶段的 SPRITE 数据量（Cluster Num）</center>

|       | All   | Sample    | Intra | F1    | FB (40%)  |
| ----  | ----  | ----      | ----  | ----  | ----:     |
|Sample 1| 2,562,880 | 2,248,559 | 1,829,574 | 1,335,255 | 1,138,141 |
|Sample 2| 2,266,648 | 1,806,029 | 1,486,294 | 1,049,433 | 932,713 |

<p><br>


## 3. 绘制样本各染色体的互作热图

首先，我们将过滤后的 SPRITE 数据由 multi-ends 互作格式转化为 double-ends 格式，转化方法为排列组合。转化后的数据量参见表 2。然后，我们用 Juicer Tools Pre 和 hic2cool 将互作数据转化为 cool 格式（互作矩阵），再使用 cooler 归一化矩阵并绘制各染色体的互作热图（图 4）。

```
脚本：
Contact_Heatmap/Prepare_Data.py
Contact_Heatmap/Heatmap_Contact.py
输入：
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
输出：
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.pairs
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.hic
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/Heatmap_Contact_MS0612-5.pdf
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.pairs
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.hic
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/Heatmap_Contact_MS0612-3.odd70.pdf
```

<center>Table. 2. 不同格式下样本的互作数量</center>

|      | cluster num（multi-ends）| contact num（double-ends）|
| ---- | ----:  | ----: | 
|Sample 1| 1,138,141 | 815,239,083 |
|Sample 2|   932,713 | 430,938,706 |

<p><br>

<img src="Doc/Fig4_Contact_Heatmap.png" alt="fig 4" />
<center>Fig. 4. 样本 4 号和 6 号染色体的互作热图</center><br>


## 4. 计算 A/B compartment

我们以 cool 格式数据作为输入，使用 cooltools 的 eigdecomp.cis_eig 函数在 1Mbp 的分辨率下对染色体进行 PCA 处理，按照与 GC 含量的相似性（Pearson）排序特征向量。其中，相似性最高的特征向量被确定为 E1。随后，我们将 E1 值大于 0 和小于 0 的染色体区间分别确定为 A 和 B compartment（图 5）。

```
脚本：
AB_Compartment/Plot_ABCompartment.py
输入：
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
输出：
AB_Compartment/Lineplot_AB_Compartment_MS0612-5.pdf
AB_Compartment/Lineplot_AB_Compartment_MS0612-3.odd70.pdf
```

<img src="Doc/Fig5_AB_Compartment.png" alt="fig 5" />
<center>Fig. 5. 样本 4 号和 6 号染色体的 AB compartment</center><br>


## 5. 对 spot 进行聚类
我们以染色体区间内的互作数量（end num）为特征，对样本区的 spot 进行 k-means 聚类。

### 5.1 构建 spot-bin 互作矩阵
我们以矩阵的格式（Sample 1：3106 $\times$ 2737，Sample 2：1853 $\times$ 2737），行为 spot id，列为 bin id，统计了 spot 中每个染色体区间（bin，1 MB）的 end num。然后，我们分别统计了 spot 和 bin 的缺失率（图 6）。然后，我们过滤了缺失率较高（Top 5%）的 spot 和 bin（图 7，Sample 1：2950 $\times$ 2451，Sample 2：1760 $\times$ 2450。

```
脚本：
KMeans/Filter_SpotBin_MissRate.py
输入：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
输出：
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
<center>Fig. 6. spot-bin 互作矩阵中 spot 和 bin 的缺失率</center><br>

<img src="Doc/Fig7_Histplot_FilterMissrate.png" alt="fig 7" />
<center>Fig. 7. 过滤后（top5%）矩阵中 spot 和 bin 的缺失率</center><br>


### 5.2 归一化及填补
首先，我们截取（clip）spot 内互作数量较为极端的 bin（Top 0.5%）。然后，使用 z-score 进行归一化，其中没有互作的 bin（end num=0）视为缺失（None）。接着，使用邻近（N=2） spot 的互作信息进行填补（图 8）。最后，因为填补过程可能会改变 spot 内互作数据的分布，脚本使用 qnorm 对 spot-bin 互作矩阵进行再次的归一化。

```
脚本：
KMeans/Normalized_Imputation.py
输入：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.matrix
输出：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
```

<center><img src="Doc/Fig8_Imputation.png" width = "200" height = "240" alt="fig 8" /></center>
<center>Fig. 8. 邻近 spot（N=2）填补中心 spot</center><br>


### 5.3 PCA 降维

我们使用 PCA 降维 spot-bin 矩阵，减少矩阵的列数，提高特征维度的信息量。我们绘制了 PCA 处理后各维度的可解释方差的积累情况（图 9）。结果显示，在 dimension=200 处曲线出现了明显的拐点，说明当 dimension>200 时，维度所能解释的方差较少。所以，我们选择 PCA 处理后前 200 维度的数据作为 spot 的特征，进行后续的聚类分析，即 Sample 1：2893 $\times$ 200，Sample 2：1705 $\times$ 200。因为 PCA 不支持数据中存在 None，所以我们剔除了填补后仍然存在缺失的 spot，故此处的 spot 数量少于 3.1。

```
脚本：
KMeans/PCA.py
输入：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.matrix
输出：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.pdf
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.pdf
```


<center><img src="Doc/Fig9_PCA.png" alt="fig 9" /></center>
<center>Fig. 9. PCA 维度与可解释方差比率的关系</center><br>


### 5.4 K-means 聚类
我们以降维后的 spot-bin 矩阵为输入，使用 K-means 方法（K=30）对 spot 进行聚类。当使用不同的 seed 时，聚类结果不太稳定，我们选择类别的空间位置及形态与器官（MOSTA 数据）较为吻合的结果。经过筛选，样本 1（MS0612-5）的参数为 K=30，seed=31；样本 2（MS0612-3.odd70）的参数为 K=30，seed=11。由于结果中一些类别内的 spot 在空间上不相邻，为了方便分析，我们将其拆分为多个类。所以，最终类别数量 > 30。其中，样本 1 被划分为 35 个类（图 10）；样本 2 被划分为 32 个类（图 11）。同时，我们将类别与器官进行了关联（图 12）。其中，Sample 2 没有找到合适的参考注释，图注为作者的推测。

```
脚本：
KMeans/KMeans.py
KMeans/Scatter_KMeans_Result.py
输入：
KMeans/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
KMeans/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.FM.norm.impu.pca.matrix
输出：
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
<center>Fig. 10. Sample 1 的 K-means 聚类结果及空间分布</center><br>

<center><img src="Doc/Fig11_KMeans_MS0612-3.odd70.split.png" width = "400" height = "600" alt="fig 11" /></center>
<center>Fig. 11. Sample 2 的 K-means 聚类结果及空间分布</center><br>

<center><img src="Doc/Fig12_KMeans.png" alt="fig 12" /></center>
<center>Fig. 12. K-means 聚类结果及对应的器官</center><br>


## 6. 计算各组织不同距离的互作频率

首先，我们根据 K-means 结果提取各类别的 SPRITE 数据，并转换为 cool 文件。然后，我们利用 cooltoos 的 expected_cis 函数计算了各组织基因组层面的互作距离与互作频率的关系，其中肝脏较为特殊（图 13）。结果显示，在不同的距离下，各组织互作频率的规律也不同。我们将互作分为短距离（≤ 10 MB）、中距离（10 MB ＜ x ≤ 40 MB）和长距离（＞ 40 MB）。短距离互作中，腹部的互作频率较高，面部互作频率中等，脑部互作频率较低；中距离互作中，肝脏互作频率较高，其他组织互作频率较低；远距离互作中，脑部及面部互作频率较高，腹部互作频率中等，肝脏互作频率较低。

```
脚本：
Contact_Distance/Prepare_TissueData.py
Contact_Distance/Prepare_CoolData.py
Contact_Distance/Lineplot_ContactDistance.py
Contact_Distance/Heatmap_ContactDistance.py
输入：
Data/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.sprite
Data/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.sprite
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
输出：
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
<center>Fig. 13. 各组织不同距离下的互作频率</center><br>


## 7. 分析各组织的区室化强度

### 7.1 计算各组织的 A/B compartment

首先，我们根据 K-means 结果提取各类别的 SPRITE 数据。然后，按照之前的步骤计算各类别各染色体的 E1 值及 A/B compartment。接着，我们分染色体计算类别的 E1 值与胚胎的 E1 值之间的 Pearson 相关系数（图 13）。结果显示一些染色体的 E1 值在类别和胚胎之间差异较大，可能是由于实验误差所致。因此，我们剔除了差异较大的染色体，使用相关性较高的 10 条染色体（2、4、5、6、8、9、11、15、17、19）进行后续分析。

```
脚本：
Compartmentalization_Strength/Prepare_CoolData.py
Compartmentalization_Strength/Plot_ABCompartment.py
输入：
Data/Tissue_Pairs/*
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
输出：
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
<center>Fig. 14. 各组织各染色体的 E1 值与胚胎 E1 值之间的相关系数</center><br>


### 7.2 计算各组织的区室化强度

我们使用 cooltools 的 saddle() 函数，将各组织的染色体按照胚胎的 E1 值划分为 40 个 bin，并计算各 bin 之间的互作频率与期望之间比值，结果以马鞍图的方式展示（图 14A）。图中 bin 在 x 轴和 y 轴上按照 E1 值的大小排序，四个角分别表示 BB、BA、AB 和 AA compartment 之间的互作频率比值。然后，我们使用 cooltools 的 saddle_strength() 函数计算各组织的区室化强度（公式如下，图 14B），其中 AA、BB、AB 和 BA 表示鞍点图四个角 5*5（extend=5）区域的互作频率比值的均值（图 14A）。结果显示，胚胎腹部的区室化强度较高，面部的区室化强度中等，脑部的区室化强度较低。

$$Strength=\frac{𝐵𝐵+AA}{𝐵A+A𝐵}$$

```
脚本：
Compartmentalization_Strength/Compartment_Strength.py
输入：
Contact_Heatmap/clusters_MS0612-5.sample.intra.FRI.F1.FB_0.4.cool
Contact_Heatmap/clusters_MS0612-3.odd70.sample.intra.FRI.F1.FB_0.4.cool
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.value
KMeans/KMeans_MS0612-5.split.label
KMeans/KMeans_MS0612-3.odd70.split.label
Compartmentalization_Strength/Data_Cool/*
输出：
Compartmentalization_Strength/Saddle.Strength_MS0612-5.pdf
Compartmentalization_Strength/Saddle.Strength_MS0612-3.odd70.pdf
Compartmentalization_Strength/Strength_MS0612-5.tsv
Compartmentalization_Strength/Strength_MS0612-3.odd70.tsv
Compartmentalization_Strength/Scatter.Strength_MS0612-5.pdf
Compartmentalization_Strength/Scatter.Strength_MS0612-3.odd70.pdf
```


<center><img src="Doc/Fig15_Compartmentalization_Strength.png" alt="fig 15" /></center>
<center>Fig. 15. 马鞍图及各组织的区室化强度</center><br>


## 8. 分析不同区室化强度组织的互作频率特征

我们分别从高、中、低区室化强度的区域内分别筛选了 sample 1 的 2、3、5 个组织（图 16A）。然后，我们以 compartment 为单位，计算了不同强度的各组织的 compartment 内和同类型 compartment 间的互作频率，并使用 bootstrapping 检验。结果显示，随着区室化强度的增加，compartment 间互作强度降低（图 16BC），compartment 内互作强度增加（图 16DE）。

```
脚本：
Compartment_Contact/Calculate_Compartment_Border.py
Compartment_Contact/Calculate_Compartment_ContactFrequency.py
Compartment_Contact/Violinplot_Compartment_ContactFrequency.py
输入：
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
Compartmentalization_Strength/E1_Value/MS0612-3.odd70.e1.value
Data/Tissue_Pairs/*
输出：
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
<center>Fig. 16. 不同区室化强度组织的互作频率分布</center><br>


## 9. 分析肝组织长距离、区室间互作的特征

因为肝组织（Label ID=1）的远距离互作（>40MB）显著低于其他组织（图 13A），所以我们具体分析了 Sample 1 的上述 10 个组织及肝组织的同类型 compartment 之间的互作差异（图 17）。图中互作颜色与互作频率正相关，互作频率越高，线条颜色越深。结果显示，肝组织的 B-B compartment 间互作频率在 2、5、9 号染色体上显著降低。

```
脚本：
Long_Compartment_Contact/Prepare_LongContact.py
Long_Compartment_Contact/Lineplot_LongContact.py
输入：
Data/Tissue_Pairs/*
Data/Tissue_Chrom_Pairs/*
Compartment_Contact/Compartment_AB_border.MS0612-5.tsv
Long_Compartment_Contact/mm10_chrom_sizes.txt
输出：
Data/Tissue_Chrom_Pairs_LongDistance/*
Long_Compartment_Contact/Compartment_ContactNum_Inter.LongDistance.MS0612-5.tsv
Long_Compartment_Contact/Compartment_ContactNum_Inter.LongDistance.MS0612-5.pdf
```


<center><img src="Doc/Fig17_LongDistance_Inter-Compartment_Contact.png" alt="fig 17" /></center>
<center>Fig. 17. 不同组织长距离、区室间互作频率的差异</center><br>


## 10. 计算各组织染色体的三维结构

我们使用 Hickit 软件计算了 11 个组织 10 条染色体的三维结构，分辨率为 1 MB（图 18）。同时，我们标注了染色体片段所属的 compartment 类别（A/B）。

```
脚本：
3D_Genome/Prepare_Data.py
3D_Genome/Scatter_3DGenome.py
3D_Genome/Scatter_3DGenome_AB.py
输入：
Data/Tissue_Chrom_Pairs/*
Compartmentalization_Strength/E1_Value/MS0612-5.e1.value
输出：
Data/Tissue_Chrom_Pairs_3DGenome/*
3D_Genome/Data_Hickit/*
3D_Genome/Scatter_3D_Structure/*
3D_Genome/Scatter_3D_Structure_AB/*
```

<center><img src="Doc/Fig18_Liver_Chrom8_3DStructure.png" alt="fig 18" /></center>
<center>Fig. 18. 肝组织（Label ID=1）8 号染色体的 3D 结构</center><br>