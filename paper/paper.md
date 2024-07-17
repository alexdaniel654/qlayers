---
title: '3DQLayers: Volumetric Layer Based Analysis for Quantitative Renal MRI'
tags:
    - Python
    - Medical Imaging
    - MRI
    - Quantitative Imaging
    - Image Analysis
authors:
    - name: Alexander J Daniel
      orcid: 0000-0003-2353-3283
      affiliation: 1
      corresponding: true
    - name: Susan T Francis
      orcid: 0000-0003-0903-7507
      affiliation: "1, 2"

affiliations: 
    - name: Sir Peter Mansfield Imaging Centre, University of Nottingham, Nottingham, United Kingdom
      index: 1
    - name: NIHR Nottingham Biomedical Research Centre, Nottingham University Hospitals NHS Trust and the University of Nottingham, Nottingham, United Kingdom
      index: 2
date: 12 July 2024
bibliography: paper.bib
---

# Summary
### Why renal MRI?
The kidneys are a pair of structurally and functionally complex organs in the lower abdomen that participate in the control of bodily fluids by regulating the balance of electrolytes, excreting waste products of metabolism and excess water from blood to urine [@lote_principles_2012]. The each kidney is separated into two tissue types; cortical tissue located towards the outside of each organ, and medullary tissue arranged in small pyramids towards the centre of the organ [@hall_guyton_2015], \autoref{fig:renal_structure}. Magnetic Resonance Imaging (MRI) is the ideal medical imaging modality to study the kidneys due to its non-ionising, non-invasive and quantitative nature [@francis_magnetic_2023; @selby_assessment_2024]. Quantitative MRI is the process of taking measurements where the value of each voxel has numerical significance, in physical units, based on the tissues underlying properties rather than simply representing signal intensity in arbitrary units. To help interpret quantitative images, regions of interest (ROI) are defined and statistical measures taken of the voxels within each region. 

### Why layers?
Segmenting ROI for the renal cortex and medulla manually is time consuming, difficult and prone to intra- and inter-reader variation thus decreasing the repeatability of measurements. Pruijm _et al_ proposed an alternative to tissue ROI based analysis in the Twelve Layer Concentric Object (TLCO) method [@piskunowicz_new_2015; @milani_reduction_2017; @li_renal_2020] where users delineate the inner and outer boundaries of the kidney to generate twelve equidistant layers between the renal pelvis and the surface of the kidney. The outer layers are analogous to the cortex and inner layers are analogous the medulla.

### Why 3DQLayers?
TLCO requires the MR image to be a single slice cutting through the kidney on their longest axis (coronal-oblique) however, this is not always desirable [@bane_consensus-based_2020]. Often researchers prefer to acquire multi-slice images to increase the number of voxels in the image and gain a better understanding of the heterogeneity of the tissue. Additionally flexibility in the orientation images are acquired at is highly desirable. These limitations of TLCO were the motivation for the development of 3DQLayers, a volumetric, quantitative-depth based analysis method for renal MRI data.
 
# Statement of need


# Figures
![a. A schematic of the kidneys showing the renal cortex and medullary pyramids. b. An anatomical MR Image of the abdomen showing the kidneys with the renal cortex appearing as a light band towards the outside of the kidneys and medullary pyramids as darker patches towards the centre of the kidneys. \label{fig:renal_structure}](kidney_overview.png){ width=90% }

# Acknowledgements

# References