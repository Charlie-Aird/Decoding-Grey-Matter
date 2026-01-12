# Decoding Grey Matter

**Decoding Grey Matter** is a framework for large-scale analysis of brain cell morphometry aimed at informing microstructural modeling of diffusion MRI (dMRI) signals.  
The project focuses on extracting detailed morphological characteristics from neural reconstructions represented as connected graphs in **SWC format**.

This repository provides MATLAB tools for analyzing neuronal structure, shape, and topology from digital reconstructions, enabling quantitative studies of neural microarchitecture.

---

## Project Overview

Neural cell morphology plays a critical role in understanding brain microstructure and its influence on diffusion MRI signals. This framework enables:

- Automated analysis of reconstructed neural cells  
- Extraction of biologically meaningful morphometric features  
- Support for large-scale datasets of neuronal reconstructions  

Neurons are represented as connected graphs and processed directly from standard **SWC** morphology files.

---

## Features

The framework extracts three main classes of morphometric descriptors:

### 1. Structural Characteristics

Fundamental biological features, including:

- Soma radius  
- Total dendritic length  
- Branch length distributions  
- Dendritic tortuosity  
- Number of branches and bifurcations  

### 2. Shape Characteristics

Orientaional characteristics of cellular processes:

- Ceullular fractional anisotropy  
- Dentrite orientaion characterised by Watson distribution   

### 3. Topological Characteristics

Topological characteristics of cellular processes via persistence diagram:

- Barcodes representing decomposition along path length
- Comparison between persistence images

## Input Format

The framework operates on neuronal reconstructions stored in **SWC format**, a standard file format for representing neuronal morphology as connected trees.
There are many open repositories hosting diverse array of cellular reconstructions such as Neuromorpho.org. 

Any valid SWC reconstruction can be analyzed using this framework.


## Requirements

This project is implemented in **MATLAB** and depends on the following:

- **TREES Toolbox for MATLAB**  
  https://www.treestoolbox.org/

- **Blender**
  https://www.blender.org/download/releases/2-79/

Make sure the toolbox and Blender are installed and added to your MATLAB path before running the code.

---

## Getting Started

Clone the repository:

```bash
git clone https://github.com/Charlie-Aird/Decoding-Grey-Matter.git
```
---

## Work flow

1. The main script, analyse_morphology.m, first identifies SWC files from the folder Reconstructions given the species and cell types being analysed.
For demonstration purposes he folder contains 10 SWC files of pyramidal cells for:
- Human [1]  
- Monkey [2]  
- Mouse [3]  
- Rat [4]  

All reconstructions were obtained from the open repository Neuromorpho.org (https://neuromorpho.org/).

Once the SWC files are identified the script has four main steps:

-**Clean trees**
 This step reads the SWC files present in the given folder and assess their integrity and removes any axonal component present within the reconstruction (this was opted for   due to the inconsistent pressence and quality of axons within the dataset)

-**Analyse Soma**
 This step quantifies metrics assocciated with the soma, such as volume and surface area, of each cell.

-**Analyse Dendrites**
 This step quantifies metrics assocciated with the dendrites, such brnach lngth, radius, and tortuosity.

-**Analyse Shape**
 This step assess the angular distribution of the dentrites around their princible axis, returning the fractional anisotropy and Watson parameter.

-**Analyse Topology**
 This step assess the dendrite topology through tropological persistence along path length from the soma. 

 2. Statistics are collated with the scripts make_stats.m and make_stats_shape.mat and saved to the folder complete_summary_staistics and can be plotted using the script StatsPlot.m found in the Images folder.

 3. Topological analysis is computed by creating persistence images for each species and cell type from the corresponding barcodes as well as visualised using the script DistancePlot.m found in the Images plot.  
    
## References

[1] Human pyramidal cells — Neuron (2016)  
https://doi.org/10.1016/j.neuron.2016.10.019  

[2] Monkey pyramidal cells — Neuron (2016)  
https://doi.org/10.1016/j.neuron.2016.02.038  

[3] Mouse pyramidal cells — Journal of Neuroinflammation (2020)  
https://doi.org/10.1186/s12974-020-01790-9  

[4] Rat pyramidal cells — Journal of Comparative Neurology  
https://doi.org/10.1002/cne.903620103  
