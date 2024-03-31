## GeoMX introduction.
GeoMx Digital Spatial Profiler (DSP) from NanoString is a state-of-the-art platform designed to address the complexities of tissue heterogeneity and the microenvironment in spatial biology. This technology is particularly advantageous for its biology-driven profiling, allowing researchers to select and analyze regions of interest with unprecedented precision and sensitivity. It uniquely enables non-destructive RNA and protein profiling from specific tissue compartments and cell populations, leveraging automated and scalable workflows compatible with FFPE and fresh frozen sections. GeoMx DSP supports various applications, from biomarker discovery to elucidating disease mechanisms and drug actions, by spatially profiling the whole transcriptome and over 570 proteins from intact tissue. It integrates with standard histology staining procedures, ensuring seamless adoption into existing research protocols. For more detailed information, visit the official NanoString website [NanoString](https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-overview/).
## Preprocessing
- **Loading Data**: This step involves the initial setup including the following key files:
    - **DCCs files**: expression count data and sequencing quality metadata
    - **PKCs file(s)**: probe assay metadata describing the gene targets present in the data
    - **Annotation file**: useful tissue information, including the type of segment profiled, segment area/nuclei count, and other tissue
- **Quality Control (QC) & Preprocessing**: This phase is crucial for ensuring data reliability. Regions of Interest (ROI)/Areas of Interest (AOI) segments and genes are selected based on quality control (QC) or limit of quantification (LOQ) metrics. The following five steps should be followed in this step.
    -	**Segment QC**: This process assesses sequencing quality and adequate tissue sampling for every ROI/AOI segment based on QC parameters like Raw sequencing reads, percentage of aligned, percentage of trimmed, etc. Those segments with low quality should be removed.
    -	**Probe QC**: This process evaluates probe performance and removes those that perform poorly. It ensures that only reliable probes contribute to gene-level count data, which is essential for accurate downstream analysis.
    -	**Create Gene-level Count Data**: This step constructs a gene-level count matrix. The count for any gene with multiple probes per segment is calculated as the geometric mean of those probes.
    -	**Limit of Quantification (LOQ)**: Determining the LOQ for each segment allows for identifying genes expressed above background levels. This step is foundational for focusing on biologically relevant data.
    -	**Filtering**: After establishing the LOQ, segments, and genes with a low signal are filtered out. This step refines the dataset further, ensuring a focus on significant, biologically relevant information.
- **Normalization**: Normalization is applied to mitigate technical variations and prepare the data for analysis. The two common methods for normalization of DSP-NGS RNA data are quartile 3 (Q3) and background normalization, which should be selected given the number of negative probe counts. If the number of negative probe counts are low, quartile 3 (3Q) is recommended. For more detailed information, please visit 
[tutorial](https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)

## GeoXM batch effect
While robust, GeoMx Digital Spatial Profiler (DSP) data often includes technical variability that can overshadow biological insights if not adequately accounted for. Batch effects, stemming from technology differences or sample preparation variations, are a significant source of such variability. These effects can lead to incorrect interpretations of biological differences, risking false discoveries. The primary goal of GeoMX datasets Batch Effect correction is to ensure the integrity of biological insights by effectively removing technical noise, enabling accurate and reliable downstream analyses.

## Batch effect solution
 - ## Tool name (StandR)
   The standR package is designed as a comprehensive R/Bioconductor tool for analyzing NanoString GeoMx Digital Spatial Profiler data, addressing limitations in current bioinformatics workflows by enhancing quality control and enabling deep insights into complex gene expression within tissues' spatial contexts. [StandR](https://academic.oup.com/nar/article/52/1/e2/7416375)
 - ## Tool principle
   The batch effect correction method in StandR is RUV4, a method for correcting unwanted variation in high-dimensional data with negative controls, operates through a series of steps aimed at distinguishing and adjusting for variations that obscure the true signal. Firstly, it utilizes negative controls—known variables unaffected by the experimental conditions—to identify unwanted variation. This unwanted variation is then mathematically separated from the data, ensuring that subsequent analyses focus on the genuine biological differences. The approach is designed to enhance the accuracy of data interpretation by effectively mitigating the impact of confounding factors, such as batch effects, without the need for precise knowledge of the number of these unwanted factors, simplifying its application across different datasets.
   For more detailed information, please visit [StandR tutorial](https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html)
 - ## Why need Grid Search
   Given the default settings, employing Batch Effect Correction (RUV4) within StandR does not efficiently discard appropriate features to distinctly separate ROI/AOI with varying cell types, while concurrently blending data from different batches to mitigate batch effects. Therefore, while retaining the framework of the StandR tool, we aim to optimize the outcomes of batch effect elimination by experimenting with various parameter configurations. In the program, three primary parameters are subject to adjustment:

    - **normalization_methods**: We intend to explore different normalization methods, including Total Count per Million (CPM), Trimmed Mean of M-values (TMM), upper quartile, and size factor.

    - **findNCG_topn**: This parameter determines the number of negative control genes selected from the top, with the adjustment range spanning from 100 to 1500 in increments of 200.

    - **BatchCorrection_k**: This refers to the count of unwanted factors to be utilized. Following the guidance in RUV's documentation, it is advisable to employ the smallest possible k where the technical variation is no longer apparent, with the parameter adjustment range being from 3 to 11 in steps of 2.

    Through conducting a grid search on the aforementioned combinations of parameters, we anticipate identifying the optimal result for batch effect correction achievable with StandR.

 - ## Evaluation
   To enhance the accuracy of downstream analyses, a successful batch effect method should preserve biological variations, such as cell types, functions, and pathological conditions, by removing unnecessary features, while eliminating influences irrelevant to the study, like differences in experimental techniques and batches. Evaluation methods include graphical representations like UMAP and t-SNE, as well as indices for clustering quality such as kBET, LISI, and ASW. In this program, we will employ the following methods to assess the effectiveness of batch effect correction.
    ### KBET
   The k-Nearest neighbor batch-effect test (kBET) is primarily utilized to assess the local mixing effect between different batches after batch correction, based on dimensionality reduction via singular value decomposition. A low rejection rate indicates effective mixing of batches (under a null hypothesis assuming uniform blending across batches, where local and overall distributions are sufficiently similar to avoid the claim of poor batch mixing). The rejection rate ranges between 0 and 1, with rates close to 0 indicating successful mitigation of batch effects, allowing cells from different batches to mix well.
   For more detailed information, please visit [KBET](https://www.nature.com/articles/s41592-018-0254-1)
    ### Silouettes
   The Silhouettes Coefficient is a metric for assessing the quality of clustering. It evaluates clustering performance by combining cohesion and separation, assigning each sample a score between 0 and 1. The overall clustering effectiveness is determined by averaging these scores across all samples. This method facilitates the comparison of different algorithms or configurations on the same dataset, examining their impact on clustering outcomes.For more detailed information, please visit
   [Silouettes](https://www.sciencedirect.com/science/article/pii/0377042787901257)

## Code demonstration
 Example code (Rmarkdown file) in code/GeoMX tutorial.Rmd
 
