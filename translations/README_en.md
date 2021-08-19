**Author**: Alessandro Giudice    
**Contributor**: Samuel Santhosh Gomez  

# Percentage of non-functionality of a patch in a unit disk
Below I report the procedure to calculate the percentage of non-functionality of a specific area of ​​a protein surface in 3D.  
The `text` written in this way represents the variables of the code used, visible in the appendix.   

## Selecting a patch
The entire protein surface studied is visible in Figure 0.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/entire_protein.png" width=600px></p>
<p align="center"><i>Figure 0</i>: The entire protein surface in 3D.</p>

A patch, like the one in Figure 1, is a group of points on a 3D surface. These points are selected as patch points if they have a mutual distance not exceeding the threshold value `Dpp`, and if they are contained in a sphere having as its center the index of a point on the surface `center` and a radius `Rs`.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point5000.png" width=450px><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point19841.png" width=450px></p>
<p align="center"><i>Figure 1</i>: Patch of point 5000 (left) and patch of point 19841 (right) of the surface.</p>

Then the selected patch must be incorporated in a cone as in Figure 2. This cone is placed along the z axis, with origin at the point C = (0,0,`z`), so that the maximum angle between the perpendicular axis and the secant connecting C to a point on the surface (or patch) is equal to `theta_max = 45`.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Cone_Point5000.png" width=600px></p>
<p align="center"><i>Figure 2</i>: The patch (red) of the point 5000 inside the cone (blue).</p>

## Creation of the fit plan
Each point of the patch is projected onto a 2D square grid with an `Npixel` side, in which each cell is a pixel. Inside each pixel there is the value of the mean or variance of the distances between the relative points of the patch and the origin C of the cone. Consequently, there are two possible grids to create:

* the matrix of the means of the distances in each pixel, as in the left part of Figures 3-4.

* the matrix of the variances of the distances in each pixel, as in the right part of Figures 3-4.

The distances used for each matrix in Figures 3-4 are only those contained in a unit disk (distances less than or equal to one).  
There are two methods to create the plane of the mean and variance of these distances.  
The first method (`CreatePlane_Weights` function) builds a grid in which the value (mean or variance) of each pixel is based on the distances between the patch points and the point C. The relative distance to a patch point ends in one pixel if the point is orthogonally above that pixel. Two examples of this method are shown in Figure 3. The top and bottom parts of the figure refer to the patch with `center = 5000` and the patch with` center = 19841` respectively. This method can be called the *Weights* method.

<p align="center">
<img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Point_5000_Weights.png" width=700px>
<img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Point_19841_Weights.png" width=700px>
</p>
<p align="center"><i>Figure 3</i>: Mean and variance of two patches (one per row) produced with the first method (<i>Weights</i> method).</p>

In the second method (`CreatePlane_Projections` function) the grid is constructed so that each pixel has a value (mean or variance) based on the distances between the patch points and the point C. Unlike the first method, the distance relative to a patch point ends in a pixel if the segment joining a patch point and point C intercepts that pixel. The same examples of Figure 3 produced with this method are shown in Figure 4. This method can be called the *Projections* method.

<p align="center">
<img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Point_5000_Projections.png" width=700px>
<img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Point_19841_Projections.png" width=700px>
</p>
<p align="center"><i>Figure 4</i>: Mean and variance of two patches (one per row) produced with the second method (<i>Projections</i> method).</p>

## Percentage of non-functionality
The percentage of non-functionality `perc` of a patch coincides with the percentage of pixels in the matrix that contain a variance greater than a `threshold`. The found value of `perc` and the value chosen for `threshold` are shown in the title of the right part of the graphs in Figures 3-4. Furthermore, the threshold value is also indicated on the relative colored bar of these figures. For each pixel, if the variance is less than this threshold, a uniform color is shown (patch with `center = 5000` in Figure 3-4), otherwise a more or less dark color is displayed for a high or low value of the variance (patch with `center = 19841` in Figure 3-4).  
The values of `perc` are calculated with the functions `PercHigherVariance_Weights` and `PercHigherVariance_Projections`. These values for each point of the surface are visible in Figure 4 and Figure 5 for the first and second method respectively.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/all_perc.png" width=800px></p>
<p align="center"><i>Figura 4</i>: Percentage of non-functionality with the <i>Weights</i> method for each point of the surface.</p>
<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/all_perc_projections.png" width=800px></p>
<p align="center"><i>Figura 5</i>: Percentage of non-functionality with the <i>Projections</i> method for each point of the surface.</p>

As shown in Figure 6-7, the second method produces fit plans with a generally lower percentage of non-functionality than the first method.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/hist_01.png" width=800px></p>
<p align="center"><i>Figura 6</i>: Histogram of the percentage of non-functionality with the <i>Weights</i> method.</p>
<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/hist_02.png" width=800px></p>
<p align="center"><i>Figura 7</i>: Histogram of the percentage of non-functionality with the <i>Projections</i> method.</p>

## Appendix
### Libraries and modules
The code written was executed with <a href="https://jupyterlab.readthedocs.io/en/stable/" target="_blank">JupyterLab</a> using `python 3.8`.  
The python modules used, installed via <a href="https://pip.pypa.io/en/stable/" target="_blank">pip</a> (including `jupyterlab`), are listed below.
```python
import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import pandas as pd
```
```python
from mayavi import mlab  
```
The `mayavi` module, specifically` mlab`, is needed to display 3D surfaces in a Qt window, so as to produce Figure 0-1-2.
While the basic libraries are  
```python
sys.path.append ("./ bin /")  
import ZernikeFunc as ZF  
import SurfaceFunc as SF  
```
written by <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.
### Parameters
