**Author**: Alessandro Giudice    
**Contributor**: Samuel Santhosh Gomez  

# Percentage of non-functionality of a patch in a unit disk
Below I report the procedure to calculate the percentage of non-functionality of a specific area of a protein surface in 3D.  
The `text` written in this way represents the variables of the code used, that is visible in the appendix or <a href="https://github.com/agds93/percentage_non_functionality/tree/main/code" target="_blank">here</a>.  

## Selecting a patch
The entire protein surface studied is visible in Figure 0.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/entire_protein.png" width=600px></p>
<p align="center"><i>Figure 0</i>: The entire protein surface in 3D.</p>

A patch, like the one in Figure 1, is a group of points on a 3D surface. These points are selected as patch points if they have a mutual distance not exceeding the threshold value `Dpp`, and if they are contained in a sphere having as its center the index of a point on the surface `center` and a radius `Rs`.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point5000.png" width=400px><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Patch_Point19841.png" width=400px></p>
<p align="center"><i>Figure 1</i>: Patch of point 5000 (left) and patch of point 19841 (right) of the surface.</p>

The selected patch must be rotated to be perpendicular to the xy plane, then it is incorporated in a cone as in Figure 2. This cone is placed along the z axis, with origin at the point C = (0,0,`z`), so that the maximum angle between the perpendicular axis and the secant connecting C to a point on the surface (or patch) is equal to `theta_max = 45`.

<p align="center"><img src="https://github.com/agds93/percentage_non_functionality/blob/main/img/Cone_Point5000.png" width=600px></p>
<p align="center"><i>Figure 2</i>: The rotated patch (red) of the point 5000 inside the cone (blue).</p>

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
sys.path.append("./ bin /")  
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
surf_a_ = pd.read_csv(surf_name_a)  
l_a = len(surf_a_["x"])
print("Npoints", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
```
where `surf_name_a` is the path to the *.dms* file used, available here .  
The array relative to the entire `surf_a` surface must be initialized as an object of the `Surface` class: 
```python
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
After loading the points of the complete surface, to plot the entire protein surface as in Figure 0 is used
```python
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
### Selecting a patch
A patch is built based on the parameters chosen via
```python
patch, _ = surf_a_obj.BuildPatch(point_pos = center, Dmin = Dpp)
```
To produce the graph in Figure 1, with `center = 5000`, we use:
```python
res1, c = SF.ConcatenateFigPlots([patch [:,: 3]])
SF.Plot3DPoints(res1[:, 0], res1[:, 1], res1[:, 2], c, 0.3)
```
To be usable, the patch in Figure 1 must be rotated (`patch` becomes `rot_patch`) so that it is perpendicular to the xy plane. This can be done with
```python 
rot_patch, rot_patch_nv = surf_a_obj.PatchReorientNew(patch, +1)
```
where the parameter `+ 1` (`-1`) indicates that the normal `rot_patch_nv` vector units are facing up (down).  
To find the `z` dimension of the origin C of the cone that encompasses the patch we use
```python
z = surf_a_obj.FindOrigin(rot_patch, 0)
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
def CreatePlane_Weights(label, patch, z_c, Np = 20) :
    rot_p = np.copy(patch)
    rot_p[:,2] -= z_c
    weigths = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2 + rot_p[:,2]**2)
    thetas = np.arctan2(rot_p[:,1], rot_p[:,0])
    dist_plane = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2)
    R = np.amax(dist_plane)*1.01
    if label == "mean" :
        plane = np.zeros((Np,Np))
    elif label == "var" :
        plane_var = np.zeros((Np,Np))
    else :
        plane = np.zeros((Np,Np))
        plane_var = np.zeros((Np,Np))
    rot_p[:,0] += R
    rot_p[:,1] -= R
    pos_plane = rot_p[:,:2]
    dR = 2.*R/Np
    rr_x = 0
    rr_y = 0
    for i in range(Np):
        rr_y = 0
        for j in range(Np):
            mask_x = np.logical_and(pos_plane[:,0]> rr_x, pos_plane[:,0]<= rr_x+dR)
            mask_y = np.logical_and(pos_plane[:,1]< -rr_y, pos_plane[:,1]>= -(rr_y+dR))
            mask = np.logical_and(mask_x, mask_y)
            if(len(weigths[mask]) > 0):
                if label == "mean" :
                    plane[j,i] = np.mean(weigths[mask])
                elif label == "var" :
                    plane_var[j,i] = np.var(weigths[mask])
                else :
                    plane[j,i] = np.mean(weigths[mask])
                    plane_var[j,i] = np.var(weigths[mask]) 
            rr_y += dR
        rr_x += dR
    if label == "mean" :
        return plane, weigths, dist_plane, thetas
    elif label == "var" :
        return plane_var, weigths, dist_plane, thetas
    else :
        return plane, plane_var, weigths, dist_plane, thetas
```
Instead to create the matrix of the `plane` mean and the `plane_var` variance with the second method we use the following function. The inputs are the same as the function used for the first method.
```python
def CreatePlane_Projections(label, patch, z_c, Np = 20) :
    rot_p = np.copy(patch)
    weigths = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2 + (rot_p[:,2]-z_c)**2)
    slope_angle = np.arcsin( (rot_p[:,2]-z_c) / weigths[:] )  
    shift = (rot_p[:,2]) / np.tan(slope_angle)
    dist_plane = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2) + shift
    R = np.amax(dist_plane)*1.01
    thetas = np.arctan2(rot_p[:,1], rot_p[:,0])
    rot_p[:,0] += shift[:]*np.cos(thetas[:])
    rot_p[:,1] += shift[:]*np.sin(thetas[:])
    if label == "mean" :
        plane = np.zeros((Np,Np))
    elif label == "var" :
        plane_var = np.zeros((Np,Np))
    else :
        plane = np.zeros((Np,Np))
        plane_var = np.zeros((Np,Np))
    rot_p[:,0] += R
    rot_p[:,1] -= R
    pos_plane = rot_p[:,:2]
    dR = 2.*R/Np
    rr_x = 0
    rr_y = 0
    for i in range(Np):
        rr_y = 0
        for j in range(Np):
            mask_x = np.logical_and(pos_plane[:,0]> rr_x, pos_plane[:,0]<= rr_x+dR)
            mask_y = np.logical_and(pos_plane[:,1]< -rr_y, pos_plane[:,1]>= -(rr_y+dR))
            mask = np.logical_and(mask_x, mask_y)
            if(len(weigths[mask]) > 0):
                if label == "mean" :
                    plane[j,i] = np.mean(weigths[mask])
                elif label == "var" :
                    plane_var[j,i] = np.var(weigths[mask])
                else :
                    plane[j,i] = np.mean(weigths[mask])
                    plane_var[j,i] = np.var(weigths[mask]) 
            rr_y += dR
        rr_x += dR
    if label == "mean" :
        return plane, weigths, dist_plane, thetas
    elif label == "var" :
        return plane_var, weigths, dist_plane, thetas
    else :
        return plane, plane_var, weigths, dist_plane, thetas
```
### Percentage of non-functionality
To calculate the percentage `perc` of variances higher than a certain `threshold` with the first method we use the following function. The input consists of:
* a label that determines whether to return, in addition to `perc`, the mean and/or variance matrix:
    * "mean" returns only average values
    * "var" produces variances only
    * "" or other gives mean values and variances
* the number of pixels per side of the `Npixel` grid.
* the `surf_a_obj` surface object.
* the `center` index of the surface point chosen as the patch center.
* the distance `Dpp` between the points of the patch.
* the chosen `threshold`.

```python
def PercHigherVariance_Weights(label, Npixel, surf_a_obj, center, Dpp, threshold) :
    patch, mask = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
    rot_patch, _ = surf_a_obj.PatchReorientNew(patch, 1)
    z = surf_a_obj.FindOrigin(rot_patch)
    if label == "var" :
        plane_var, _, _, _ = CreatePlane_Weigths("var", patch=rot_patch, z_c=z , Np=Npixel)
    else :
        plane, plane_var, _, _, _ = CreatePlane_Weigths("", patch=rot_patch, z_c=z , Np=Npixel)
    ZernikeM = ZF.Zernike2d(plane_var)
    plane_var_polar = ZernikeM.r
    polar_mask = (plane_var_polar <= 1)
    plane_var_masked = plane_var[polar_mask]
    Npixel_new = len(plane_var_masked)
    num_high_var = np.count_nonzero( plane_var_masked > threshold )
    perc = num_high_var / Npixel_new 
    if label == "var" :
        return plane_var, perc
    else :
        return plane, plane_var, perc
```
Instead to calculate the percentage `perc` of variances higher than a certain threshold with the second method we use the following function. The inputs are the same as the function used for the first method, in fact the only difference is the use of `CreatePlane_Projections` instead of` CreatePlane_Weights`.
```python
def PercHigherVariance_Projections(label, Npixel, surf_a_obj, center, Dpp, threshold) :
    patch, mask = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
    rot_patch, _ = surf_a_obj.PatchReorientNew(patch, 1)
    z = surf_a_obj.FindOrigin(rot_patch)
    if label == "var" :
        plane_var, _, _, _ = CreatePlane_Projections("var", patch=rot_patch, z_c=z , Np=Npixel)
    else :
        plane, plane_var, _, _, _ = CreatePlane_Projections("", patch=rot_patch, z_c=z , Np=Npixel)
    ZernikeM = ZF.Zernike2d(plane_var)
    plane_var_polar = ZernikeM.r
    polar_mask = (plane_var_polar <= 1)
    plane_var_masked = plane_var[polar_mask]
    Npixel_new = len(plane_var_masked)
    num_high_var = np.count_nonzero( plane_var_masked > threshold )
    perc = num_high_var / Npixel_new 
    if label == "var" :
        return plane_var, perc
    else :
        return plane, plane_var, perc
```
To calculate the `perc` of each point on the surface with the first method we use
```python
limit = l_a
step = 1
points_list = np.arange(0,limit,step)
perc = np.zeros((len(points_list)))
for i in range(len(points_list)) :
    _, _, perc[i] = PercHigherVariance_Weights(Npixel, Rs, surf_a_obj, points_list[i], Dpp, threshold)  
print("Number of patches =",len(perc))
```
To calculate the `percent` of each point on the surface with the second method we use
```python
limit = l_a
step = 1
points_list = np.arange(0,limit,step)
perc = np.zeros((len(points_list)))
for i in range(len(points_list)) :
    _, _, perc[i] = PercHigherVariance_Projections(Npixel, Rs, surf_a_obj, points_list[i], Dpp, threshold)
print("Number of patches =",len(perc))
```
This code takes a long time to run so the results are saved to a file:
* values found with `PercHigherVariance_Weights` are visible <a href="/data/all_perc.txt" target="_blank">here</a>.
* values found with `PercHigherVariance_Projections` are visible <a href="/data/all_perc_projections.txt" target="_blank">here</a>.
```python
with open("all_perc.txt", "w") as file0 :
    for i in range(len(points_list)) :
        file0.write("{}\t{}\n".format(points_list[i],perc[i]))
```
thus it is possible to load patch indexes and relative percentages directly from the file.
```python
points_list = np.loadtxt("./risultati/all_perc.txt", usecols=0, unpack=True)
perc = np.loadtxt("./risultati/all_perc.txt", usecols=1, unpack=True)
```
To find the highs and lows of `perc` in a specified interval, use 
```python
def zone_extremes(vect, begin, end) :
    if begin == 0 :
        border = 0
    else :
        border = begin
    if end == len(vect) :
        end -= 1
    vect_lim = vect[begin:end]
    imax = np.argmax(vect_lim) + border
    imin = np.argmin(vect_lim) + border
    vmax = np.amax(vect_lim)
    vmin = np.amin(vect_lim)
    return imax, vmax, imin, vmin
```
### Utilities
This is followed by the function that returns an array with non-null elements only in cells where the mask is True. It is needed in the `PlotMeanVariancePatch` to plot a uniform color for variance values below the chosen threshold.
```python
def MatrixMasked(matrix, mask):
    matrix_new = np.where(mask, matrix, 0)
    return matrix_new
```
### Other Charts
The following function provides the graphs of Figures 3-4. The inputs are:
* the index of the center point of the `center` patch.
* the distance between the points of the `Dpp` patch.
* the radius `Rs` of the sphere that includes the patch.
* the found value of `perc`.
* the chosen threshold value.
* the matrix of the mean.
* the variance matrix.
* the color maps to use.
* the name of the output file with the appropriate extension.
```python
def PlotMeanVariancePatch(center, Dpp, Rs, perc, T, pm, pv, color_maps, name) :
    mx1 = pm
    mx2 = MatrixMasked(pv, (pv >= T))
    matrix = [ mx1, mx2 ]
    s0 = "Patch of center = {}, distance between points = {}, radius = {}".format(center, Dpp, Rs)
    s1 = "Mean in function of position\n\n"
    if T == 0 :
        s2 = "Variance in function of position\n\n"
    else :
        s2 = "Variance in function of position.  Threshold = {:.2f}\nPercentage of higher variances = {:.4f}\n".format(T, perc)   
    titles = [s1, s2]
    if len(color_maps) != 2 :
        color_maps = ["Greens", "Reds"]
    fig, ax = mpl.subplots(nrows=1, ncols=2, figsize=(8,4), facecolor="white", dpi=200)
    fig.suptitle(s0, fontsize="9")
    for row in range(1) :
        for col in range(2):
            data = matrix[col]
            min_data = np.amin(data)
            max_data = np.amax(data)
            ax[col].set_title(titles[col], fontsize="8")
            ax[col].set_xlabel("x", fontsize="8")
            ax[col].set_ylabel("y", fontsize="8")
            ax[col].tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in ax[col].spines.keys():  # 'top', 'bottom', 'left', 'right'
                ax[col].spines[side].set_linewidth(0.30)
                ax[col].spines[side].set_color("black")
            im = ax[col].pcolormesh(data, cmap=color_maps[col], rasterized=True)
            if col == 1 :
                ticks_list = [min_data, T, max_data]
            else :
                ticks_list = [min_data, max_data]
            cb = mpl.colorbar(im, ax=ax[col], ticks=ticks_list)
            cb.ax.tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in cb.ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
                cb.ax.spines[side].set_linewidth(0.30)
                cb.ax.spines[side].set_color("black")
    fig.tight_layout()
    if name != "" or name == "default" :
        if name == "default" :
            n = "MeanVariance_Patch{}_perc{}.pdf".format(center, Dpp, perc)
        else :
            n = name
        mpl.savefig("{}".format(n))
        print("The figure was generated.")
```
The graph in Figure 4-5 is produced by
```python
fig, ax = mpl.subplots(nrows=1, ncols=1, figsize=(8,4), facecolor="white", dpi=200)
ax.set_xlim(0, len(perc))
ax.set_ylim(0, np.amax(perc)+0.01)
ax.set_title("Threshold = {}, Points = {}, Pixels = {}, Dpp = {}, Rs = {}".format(threshold,len(points_list),Npixel,Dpp,Rs), fontsize="8")
ax.set_xlabel("Surface point", fontsize="8")
ax.set_ylabel("Percentage", fontsize="8")
ax.tick_params(axis="both", width ="0.30", color="black", labelsize="6")
ax.locator_params(axis="x", nbins=21)
ax.locator_params(axis="y", nbins=21)
for side in ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax.spines[side].set_linewidth(0.30)
    ax.spines[side].set_color("black")
ax.plot(points_list, perc, "o", markersize="0.4", rasterized=True)
fig.tight_layout()
mpl.savefig("all_perc.pdf")
```
The graph in Figure 6-7 is produced by
```python
fig, ax = mpl.subplots(nrows=1, ncols=1, figsize=(8,4), facecolor="white", dpi=200)
ax.set_title("Threshold = {} $\AA$, Points = {}, Pixels = {}, Dpp = {}, Rs = {}".format(threshold,len(points_list),Npixel,Dpp,Rs), fontsize="8")
ax.set_xlabel("Percentage", fontsize="8")
ax.set_ylabel("Number of values", fontsize="8")
ax.tick_params(axis="both", width ="0.30", color="black", labelsize="6")
ax.locator_params(axis="x", nbins=20)
ax.locator_params(axis="y", nbins=20)
for side in ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax.spines[side].set_linewidth(0.30)
    ax.spines[side].set_color("black")
ax.hist(perc, bins=int(np.sqrt(len(perc))), histtype="step", rasterized=True)
fig.tight_layout()
mpl.savefig("hist.pdf")
```
