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
