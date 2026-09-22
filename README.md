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

** Usage: **

1. Install **MATLAB R2023a or later**.
2. Ensure the following MATLAB toolboxes are installed:
   - **Statistics and Machine Learning Toolbox** – used to generate the Sobol quasi-random sequence for CVT seed initialisation. If this toolbox is not available, the Sobol sequence can be replaced with standard random number generation.
   - **Signal Processing Toolbox** – used to generate the Gaussian window. If this toolbox is not available, the Gaussian window can be defined manually.
3. Add the repository and its subdirectories to the MATLAB path.
4. Open MATLAB and run:

```matlab
demo
**Citation:**

If you use this method in your research, we would greatly appreciate if you cite our original paper! The citation details can be found in the [Paper URL].

