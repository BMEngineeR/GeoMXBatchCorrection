## GeoMX introduction.
GeoMx Digital Spatial Profiler (DSP) from NanoString is a state-of-the-art platform designed to address the complexities of tissue heterogeneity and the microenvironment in spatial biology. This technology is particularly advantageous for its biology-driven profiling, allowing researchers to select and analyze regions of interest with unprecedented precision and sensitivity. It uniquely enables non-destructive RNA and protein profiling from specific tissue compartments and cell populations, leveraging automated and scalable workflows compatible with FFPE and fresh frozen sections. GeoMx DSP supports various applications, from biomarker discovery to elucidating disease mechanisms and drug actions, by spatially profiling the whole transcriptome and over 570 proteins from intact tissue. It integrates with standard histology staining procedures, ensuring seamless adoption into existing research protocols. For more detailed information, visit the official NanoString website [NanoString](https://nanostring.com/products/geomx-digital-spatial-profiler/geomx-dsp-overview/).
## Preprocessing
•	**Loading Data**: This step involves the initial setup including the following key files:
<br/> &emsp;
 •	**DCCs files** - expression count data and sequencing quality metadata
 <br/> &emsp;
 •	**PKCs file(s)** - probe assay metadata describing the gene targets present in the data
 <br/> &emsp;
 •	**Annotation file** - useful tissue information, including the type of segment profiled, segment area/nuclei count, and other tissue
 <br/>
 <br/>
•	**Quality Control (QC) & Preprocessing**: This phase is crucial for ensuring data reliability. ROI/AOI segments and genes are selected based on quality control (QC) or limit of quantification (LOQ) metrics. The following five steps should be followed in this step.
<br/> &emsp;
 •	**Segment QC**: This process assesses sequencing quality and adequate tissue sampling for every ROI/AOI segment based on QC parameters like Raw sequencing reads, percentage of aligned, percentage of trimmed, etc. Those segments with low quality should be removed.
 <br/> &emsp;
 •	**Probe QC**: This process evaluates probe performance and removes those that perform poorly. It ensures that only reliable probes contribute to gene-level count data, which is essential for accurate downstream analysis.
 <br/> &emsp;
 •	**Create Gene-level Count Data**: This step constructs a gene-level count matrix. The count for any gene with multiple probes per segment is calculated as the geometric mean of those probes.
 <br/> &emsp;
 •	**Limit of Quantification (LOQ)**: Determining the LOQ for each segment allows for identifying genes expressed above background levels. This step is foundational for focusing on biologically relevant data.
 <br/> &emsp;
•	**Filtering**: After establishing the LOQ, segments, and genes with a low signal are filtered out. This step refines the dataset further, ensuring a focus on significant, biologically relevant information.
<br/>
<br/>
•	**Normalization**: Normalization is applied to mitigate technical variations and prepare the data for analysis. The two common methods for normalization of DSP-NGS RNA data are quartile 3(Q3) and background normalization, which should be selected given the number of negative probe counts. If the number of negative probe counts are low, quartile (3Q) is recommended. 
<br/>
<br/>
For more detailed information, visit [tutorial](https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)
## GeoXM batch effect
While robust, GeoMx Digital Spatial Profiler (DSP) data often includes technical variability that can overshadow biological insights if not adequately accounted for. Batch effects, stemming from technology differences or sample preparation variations, are a significant source of such variability. These effects can lead to incorrect interpretations of biological differences, risking false discoveries. The primary goal of GeoMX datasets Batch Effect correction is to ensure the integrity of biological insights by effectively removing technical noise, enabling accurate and reliable downstream analyses.
## Batch effect solution
 - ## Tool name (StandR)
   [StandR tutorial](https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html)
 - ## Tool principle
 - ## Why need Grid Search
 - ## Evaluation
    remain Biology variation (cell type, condition, blablabla)
    reduce unwanted variation (individual differences, technology differences)
    ### KBET (paper link)
    ### Silouette (Paper link)
    
## Code demonstration
 ()[link!]
 
