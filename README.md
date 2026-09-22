# Density Encoded Line Integral Convolution (DELIC)

This code provides a powerful method for visualizing vector fields, specifically optical axis fields from polarimetry, using the Density-Encoded Line Integral Convolution (DELIC) algorithm. DELIC extends the classic Line Integral Convolution (LIC) method, enabling high-density vector field visualization with minimal clutter. The algorithm is versatile and can be applied to any general vector field for tractography.

**Features:**

- Visualization of Polarimetry Data: DELIC is designed to visualize optical axis fields, such as those from a Polarization Sensitive Optical Coherence Tomography (PS-OCT) system.
- High-Density Tractography: Visualizes complex structures like collagen fibers in skin at high densities with reduced visual clutter.
- Customizable for Other Fields: The algorithm is flexible and can be adapted to visualize other vector fields.

**Example Data:**

Example data is provided, showcasing the orientation and organization of collagen fibers in human dorsal hand skin. Data is acquired from a PS-OCT system. Additional data can be requested if needed.

**Overview:**

The demo.m script demonstrates how to use the DELIC function. This script serves as a comprehensive example, showing how to generate DELIC images from vector field data.

Key functions in this repository include:
1) DELIC.m: The main function to create DELIC images from vector fields.
2) GenerateColorDelic.m: A helper function that allows you to embed colors into your DELIC images, enhancing visualization.
3) WeightedCentroidalVoronoi.m: Computes a weighted Centroidal Voronoi Tessellation (CVT) using Lloyd's algorithm. This function can be modified for various applications beyond DELIC.

**Usage:**

1) Install MATLAB R2023a or later.
2) Ensure the Statistics and Machine Learning Toolbox and Signal Processing Toolbox are installed.
   The Statistics and Machine Learning Toolbox is required to generate a Sobel quasi-random sequence. Standard random number generation could be utilized if the toolbox is not present in your      installation
   The Signal Processing Toolbox is needed to generate a gaussian window. A gaussian window could be manually defined if the toolbox is not present in your installation.
4) Add the repository and its subdirectories to the MATLAB path.
5) Open MATLAB and run demo.m

**Citation:**

If you use this method in your research, we would greatly appreciate if you cite our original paper! The citation details can be found in the [Paper URL].

