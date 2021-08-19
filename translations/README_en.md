**Author**: Alessandro Giudice    
**Contributor**: Samuel Santhosh Gomez  

# Percentage of non-functionality of a patch in a unit disk
Below I report the procedure to calculate the percentage of non-functionality of a specific area of ​​a protein surface in 3D.  
The `text` written in this way represents the variables of the code used, that is visible in the appendix or <a href="https://github.com/agds93/percentage_non_functionality/tree/main/code" target="_blank">here</a>.  

## Selecting a patch
The entire protein surface studied is visible in Figure 0.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/entire_protein.png" width=600px></p>
<p align="center"><i>Figure 0</i>: The entire protein surface in 3D.</p>

A patch, like the one in Figure 1, is a group of points on a 3D surface. These points are selected as patch points if they have a mutual distance not exceeding the threshold value `Dpp`, and if they are contained in a sphere having as its center the index of a point on the surface `center` and a radius `Rs`.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point5000.png" width=440px><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point19841.png" width=440px></p>
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
The parameter values used to select a patch and create the fit plan are
```python
Npixel = 25    # the side of the plane in pixels
Dpp = 0.5      # the distance between points of the same patch
Rs = 6         # the radius of the sphere that includes the patch
threshold = 5  # threshold value to determine if the variance is high
```
Values are in ångström units, except `Npixel`.  
The patch center point index used for Figure 1, Figure 2, and for the graphs at the top of Figure 3-4 is
```python
center = 5000
```
instead the one used for the graphs at the bottom of Figure 3-4 is
```python
center = 19841
```
### Load the protein surface
Once the protein to be studied has been chosen, the relative *.pdb* file must be downloaded from the <a href="https://www.rcsb.org/" target="_blank">Protein Data Bank</a> from which to create a *.dms* file (for example with the tool <a href="https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html" target="_blank">dms</a>), containing a series of data on atoms and points on the surface.
```python
surf_name_a = "./data/4bs2_RRM2.dms"
surf_a_ = pd.read_csv (surf_name_a)  
l_a = len (surf_a _ ["x"])
print ("Npoints", l_a)
surf_a = np.zeros ((l_a, 6))
surf_a [:,:] = surf_a _ [["x", "y", "z", "Nx", "Ny", "Nz"]]
```
where `surf_name_a` is the path to the *.dms* file used, available here .  
The array relative to the entire `surf_a` surface must be initialized as an object of the `Surface` class: 
```python
surf_a_obj = SF.Surface (surf_a [:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
After loading the points of the complete surface, to plot the entire protein surface as in Figure 0 is used
```python
res1, c = SF.ConcatenateFigPlots ([surf_a_obj.surface [:,: 3]])
SF.Plot3DPoints (res1 [:, 0], res1 [:, 1], res1 [:, 2], c, 0.3)
```
### Selecting a patch
A patch is built based on the parameters chosen via
```python
patch, _ = surf_a_obj.BuildPatch (point_pos = center, Dmin = Dpp)
```
To produce the graph in Figure 1, with `center = 5000`, we use:
```python
res1, c = SF.ConcatenateFigPlots ([patch [:,: 3]])
SF.Plot3DPoints (res1 [:, 0], res1 [:, 1], res1 [:, 2], c, 0.3)
```
To be usable, the patch in Figure 1 must be rotated (`patch` becomes` rot_patch`) so that it is perpendicular to the xy plane. This can be done with
```python 
rot_patch, rot_patch_nv = surf_a_obj.PatchReorientNew (patch, +1)
```
where the parameter `+ 1` (` -1`) indicates that the normal `rot_patch_nv` vector units are facing up (down).  
To find the `z` dimension of the origin C of the cone that encompasses the patch we use
```python
z = surf_a_obj.FindOrigin (rot_patch, 0)
```
where replacing `0` with` 1` produces the graph of the patch embedded inside the cone. Figure 2 is obtained with `center = 5000`.
### Creation of the fit plan
To create the matrix of the `plane` mean and the `plane_var` variance with the first method, use the following function. The input consists of:
* a label that determines whether to return the mean and/or variance matrix:
    * "mean" returns only average values
    * "var" produces variances only
    * "" or other gives mean and variances values
* the patch rotated.
* the `z` dimension of point C.
* the number of pixels per side of the `Npixel` grid.
```python
def CreatePlane_Weights (label, patch, z_c, Np = 20):
    rot_p = np.copy (patch)
    rot_p [:, 2] - = z_c
    weigths = np.sqrt (rot_p [:, 0] ** 2 + rot_p [:, 1] ** 2 + rot_p [:, 2] ** 2)
    thetas = np.arctan2 (rot_p [:, 1], rot_p [:, 0])
    dist_plane = np.sqrt (rot_p [:, 0] ** 2 + rot_p [:, 1] ** 2)
    R = np.amax (dist_plane) * 1.01
    if label == "mean":
        plane = np.zeros ((Np, Np))
    elif label == "var":
        plane_var = np.zeros ((Np, Np))
    else:
        plane = np.zeros ((Np, Np))
        plane_var = np.zeros ((Np, Np))
    rot_p [:, 0] + = R
    rot_p [:, 1] - = R
    pos_plane = rot_p [:,: 2]
    dR = 2. * R / Np
    rr_x = 0
    rr_y = 0
    for i in range (Np):
        rr_y = 0
        for j in range (Np):
            mask_x = np.logical_and (pos_plane [:, 0]> rr_x, pos_plane [:, 0] <= rr_x + dR)
            mask_y = np.logical_and (pos_plane [:, 1] <-rr_y, pos_plane [:, 1]> = - (rr_y + dR))
            mask = np.logical_and (mask_x, mask_y)
            if (len (weigths [mask])> 0):
                if label == "mean":
                    plane [j, i] = np.mean (weigths [mask])
                elif label == "var":
                    plane_var [j, i] = np.var (weigths [mask])
                else:
                    plane [j, i] = np.mean (weigths [mask])
                    plane_var [j, i] = np.var (weigths [mask]) 
            rr_y + = dR
        rr_x + = dR
    if label == "mean":
        return plane, weigths, dist_plane, thetas
    elif label == "var":
        return plane_var, weigths, dist_plane, thetas
    else:
        return plane, plane_var, weigths, dist_plane, thetas
```
Instead to create the matrix of the `plane` mean and the `plane_var` variance with the second method we use the following function. The inputs are the same as the function used for the first method.
```python
 def CreatePlane_Projections (label, patch, z_c, Np = 20):
    rot_p = np.copy (patch)
    weigths = np.sqrt (rot_p [:, 0] ** 2 + rot_p [:, 1] ** 2 + (rot_p [:, 2] -z_c) ** 2)
    slope_angle = np.arcsin ((rot_p [:, 2] -z_c) / weigths [:])  
    shift = (rot_p [:, 2]) / np.tan (slope_angle)
    dist_plane = np.sqrt (rot_p [:, 0] ** 2 + rot_p [:, 1] ** 2) + shift
    R = np.amax (dist_plane) * 1.01
    thetas = np.arctan2 (rot_p [:, 1], rot_p [:, 0])
    rot_p [:, 0] + = shift [:] * np.cos (thetas [:])
    rot_p [:, 1] + = shift [:] * np.sin (thetas [:])
    if label == "mean":
        plane = np.zeros ((Np, Np))
    elif label == "var":
        plane_var = np.zeros ((Np, Np))
    else:
        plane = np.zeros ((Np, Np))
        plane_var = np.zeros ((Np, Np))
    rot_p [:, 0] + = R
    rot_p [:, 1] - = R
    pos_plane = rot_p [:,: 2]
    dR = 2. * R / Np
    rr_x = 0
    rr_y = 0
    for i in range (Np):
        rr_y = 0
        for j in range (Np):
            mask_x = np.logical_and (pos_plane [:, 0]> rr_x, pos_plane [:, 0] <= rr_x + dR)
            mask_y = np.logical_and (pos_plane [:, 1] <-rr_y, pos_plane [:, 1]> = - (rr_y + dR))
            mask = np.logical_and (mask_x, mask_y)
            if (len (weigths [mask])> 0):
                if label == "mean":
                    plane [j, i] = np.mean (weigths [mask])
                elif label == "var":
                    plane_var [j, i] = np.var (weigths [mask])
                else:
                    plane [j, i] = np.mean (weigths [mask])
                    plane_var [j, i] = np.var (weigths [mask]) 
            rr_y + = dR
        rr_x + = dR
    if label == "mean":
        return plane, weigths, dist_plane, thetas
    elif label == "var":
        return plane_var, weigths, dist_plane, thetas
    else:
        return plane, plane_var, weigths, dist_plane, thetas
```
