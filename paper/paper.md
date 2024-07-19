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
Informative measurements of the kidneys structure and function can be performed using quantitative Magnetic Resonance Imaging (MRI) where each voxel of the image is a measurement of the physical properties of the tissue being image. Traditionally, analysis of these images is performed by segmenting the kidney into its constituent tissue types and calculating the average of each measurement for each tissue type. The process of segmenting renal tissue types is time consuming and inaccurate.
An alternative to tissue segmentation proposed by Pruijm _et al_ involves dividing the kidney into layers based on each voxels distance between the outer and inner surface of the kidney, a method known as Twelve Layer Concentric Objects (TLCO) [@piskunowicz_new_2015; @milani_reduction_2017; @li_renal_2020]. Layer based analysis only requires segmenting the whole kidney rather than the tissues within and is therefore quicker and more repeatable. The TLCO method does however have some limitation, it can only be performed on a single slice image, requires the image to be acquired at a specific angle to the kidneys, requires manual delineation of the outside and inside surface of the kidney, divides the kidney into the same number of layers irrespective of the size of the kidney and the software itself is closed source.
These limitations are addressed by `3DQLayers`, a Python package to automatically define 3D, multi-slice, renal layers of known thickness.

# Statement of need
## Background
The kidneys are a pair of structurally and functionally complex organs in the lower abdomen that participate in the control of bodily fluids by regulating the balance of electrolytes, excreting waste products of metabolism and excess water from blood to urine [@lote_principles_2012]. The each kidney is separated into two tissue types; cortical tissue located towards the outside of each organ, and medullary tissue arranged in small pyramids towards the centre of the organ [@hall_guyton_2015], as shown in \autoref{fig:renal_structure}. Quantitative MRI is the process of taking measurements where the value of each voxel has numerical significance, in physical units, based on the tissues underlying properties rather than simply representing signal intensity in arbitrary units. Example quantitative measurements include how readily water can diffuse through the tissue and the rate at which blood perfuses into the tissue. To help interpret quantitative images, regions of interest (ROI) are defined and statistical measures taken of the voxels within each region. 

Segmenting ROI for the renal cortex and medulla manually is time consuming, difficult and prone to intra- and inter-reader variation thus decreasing the repeatability of measurements. Pruijm _et al_ proposed an alternative to tissue ROI based analysis in the Twelve Layer Concentric Object (TLCO) method [@piskunowicz_new_2015; @milani_reduction_2017; @li_renal_2020] where users delineate the inner and outer boundaries of the kidney to generate twelve equidistant layers between the renal pelvis and the surface of the kidney. The outer layers are analogous to the cortex; inner layers, the medulla; and gradient of the central layers, the cortico-medullary difference.

TLCO requires the MR image to be a single slice cutting through the kidneys on their longest axis (coronal-oblique) however, this is not always desirable [@bane_consensus-based_2020]. Often researchers prefer to acquire multi-slice images to increase the number of voxels in the image and gain a better understanding of the heterogeneity of the kidney. Additionally flexibility in the orientation images are acquired at is highly desirable. These limitations of TLCO were the motivation for the development of `3DQLayers`, a volumetric, quantitative-depth based analysis method for renal MRI data.

## Methods
`3DQLayers` is an open-source Python package that aims to build upon the premise of TLCO and allow layer based analysis to be fully automated for use in large studies. `3DQLayers` fundamentally differs from TLCO in that layers are defined based on each voxels distance from the surface of the kidney in millimetres rather than the proportion through the kidney. As such, the input to `3DQLayers` is a whole kidney ROI, this can be automatically generated from a structural image [@daniel_automated_2021; @daniel_renal_2024].
The pipeline by which layers are defined is outlined in \autoref{fig:flowchart}. Pre-processing steps fill in holes in the ROI caused by cysts as the surface of a cyst is not representative of the surface of the kidney. Next the voxel-based representation of the ROI is converted to a smoothed mesh-based representation of the kidneys, the distance from the centre of each voxel to the surface of the mesh can then be calculated producing a depth map [@dawson-haggerty_trimesh_2023]. As tissue adjacent to the renal pelvis is not representative of the medulla, it is excluded from layer-based analysis. This is achieved by automatically segmenting the pelvis then calculating the distance from each voxel to the pelvis as above. Voxels closer than a specified threshold, typically 10 mm, are excluded from the depth map. Finally, a layer image is generated by quantising the depth map to a desired layer thickness, typically 1 mm. 
The layer image and quantitative images are resampled to the same resolution using `NiBabel` [@brett_nipy/nibabel_2019], this allows each layer to be used as an ROI with statistical measures of the quantitative image e.g. median, standard deviation and kurtosis, calculated as a function of depth through the kidney. The gradient of the central layers can be calculated, additionally, if tissue ROI are available the distribution of tissue types with depth can be explored.
An object oriented interface makes it easy for users to generate layers and use them to analyse quantitative images. [Documentation](https://qlayers.readthedocs.io/) is provided to guide users through instillation via `PyPI`, `conda` or from [source code on GitHub](https://github.com/alexdaniel654/qlayers) and includes tutorials and an API reference. An automated test suite with high coverage gives users confidence in the stability of `3DQLayers` and that there will be no unexpected changes to results unless highlighted in the change-log.
## Usage Examples


# Figures
![a) A schematic of the kidneys showing the renal cortex and medullary pyramids. b) An anatomical MR Image of the abdomen showing the kidneys with the renal cortex appearing as a light band towards the outside of the kidneys and medullary pyramids as darker patches towards the centre of the kidneys. \label{fig:renal_structure}](kidney_overview.png){ width=90% }

![The mask from the T2-weighted structural scan (a i) has any cysts filled (a ii) and is converted into a smooth mesh representing the renal surface (b i and ii). The distance (in mm) from each voxel to the surface of the mesh is calculated (b iii). The renal pelvis is segmented (c i) and any tissue within 10 mm (c ii) of the pelvis is excluded from the depth map (c iii). The tissue is then grouped into layers of a desired thickness, here shown as 5 mm layers for illustrative purposes (d). \label{fig:flowchart}](flowchart.png)

# Acknowledgements

# References