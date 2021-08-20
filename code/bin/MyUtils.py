import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import pandas as pd
from mayavi import mlab

sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF


# Modified version of the CreatePlane method (see SurfaceFunc.py line 893).
# it, in addition to the plan with the average values, also returns the plan with the variances of the chosen patch.
def CreatePlane_Weights(label, patch, z_c, Np = 20) :

    rot_p = np.copy(patch)

    # shifting patch to have the cone origin in [0,0,0]
    rot_p[:,2] -= z_c

    # computing distances between points and the origin
    Weights = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2 + rot_p[:,2]**2)
    
    # computing angles on plane
    thetas = np.arctan2(rot_p[:,1], rot_p[:,0])

    # computing distances in plane
    dist_plane = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2)

    # computing the circle radius as the maximum distant point
    R = np.amax(dist_plane)*1.01

    # creating plane matrix
    if label == "mean" :
        plane = np.zeros((Np,Np))
    elif label == "var" :
        plane_var = np.zeros((Np,Np))
    else :
        plane = np.zeros((Np,Np))
        plane_var = np.zeros((Np,Np))
    
    # adapting points to pixels
    rot_p[:,0] += R
    rot_p[:,1] -= R

    pos_plane = rot_p[:,:2]
    
    # taking values only within unit disk

    dR = 2.*R/Np
    rr_x = 0
    rr_y = 0
    for i in range(Np):
        rr_y = 0
        for j in range(Np):
            mask_x = np.logical_and(pos_plane[:,0]> rr_x, pos_plane[:,0]<= rr_x+dR)
            mask_y = np.logical_and(pos_plane[:,1]< -rr_y, pos_plane[:,1]>= -(rr_y+dR))
            mask = np.logical_and(mask_x, mask_y)
            if(len(Weights[mask]) > 0):
                if label == "mean" :
                    plane[j,i] = np.mean(Weights[mask])
                elif label == "var" :
                    plane_var[j,i] = np.var(Weights[mask])
                else :
                    plane[j,i] = np.mean(Weights[mask])
                    plane_var[j,i] = np.var(Weights[mask]) 
            rr_y += dR
        rr_x += dR
        
    if label == "mean" :
        return plane, Weights, dist_plane, thetas
    elif label == "var" :
        return plane_var, Weights, dist_plane, thetas
    else :
        return plane, plane_var, Weights, dist_plane, thetas




# Alternate version of CreatePlane_Weights.
# This time the distances are collected based on the intersections of the segments with that plane.
# Each segment starts at a point on the surface and ends at the cone origin. 
def CreatePlane_Projections(label, patch, z_c, Np = 20) :
    
    rot_p = np.copy(patch)
     
    # computing distances between points and the origin of cone in [0,0,z_c]
    Weights = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2 + (rot_p[:,2]-z_c)**2)
    
    # angles of slope, respect to the fit plane, of segments that connect a surface point and the origin of cone
    slope_angle = np.arcsin( (rot_p[:,2]-z_c) / Weights[:] )
    
    # shift from original xy coordinates to intersection coordinates of segments with fit plane
    # shift values are positive if rot_p[:,2] > 0 or negative if rot_p[:,2] < 0  
    shift = (rot_p[:,2]) / np.tan(slope_angle)
    
    # computing distances in plane
    dist_plane = np.sqrt(rot_p[:,0]**2 + rot_p[:,1]**2) + shift
    
    # computing the circle radius as the maximum distant point
    R = np.amax(dist_plane)*1.01
    
    # computing angles on plane (slope of shift segments on fit plane)
    thetas = np.arctan2(rot_p[:,1], rot_p[:,0])
    
    # new coordinates of points on fit plane
    rot_p[:,0] += shift[:]*np.cos(thetas[:])
    rot_p[:,1] += shift[:]*np.sin(thetas[:])

    # creating plane matrix
    if label == "mean" :
        plane = np.zeros((Np,Np))
    elif label == "var" :
        plane_var = np.zeros((Np,Np))
    else :
        plane = np.zeros((Np,Np))
        plane_var = np.zeros((Np,Np))
    
    # adapting points to pixels
    rot_p[:,0] += R
    rot_p[:,1] -= R

    pos_plane = rot_p[:,:2]
    
    # taking values only within unit disk

    dR = 2.*R/Np
    rr_x = 0
    rr_y = 0
    for i in range(Np):
        rr_y = 0
        for j in range(Np):
            mask_x = np.logical_and(pos_plane[:,0]> rr_x, pos_plane[:,0]<= rr_x+dR)
            mask_y = np.logical_and(pos_plane[:,1]< -rr_y, pos_plane[:,1]>= -(rr_y+dR))
            mask = np.logical_and(mask_x, mask_y)
            if(len(Weights[mask]) > 0):
                if label == "mean" :
                    plane[j,i] = np.mean(Weights[mask])
                elif label == "var" :
                    plane_var[j,i] = np.var(Weights[mask])
                else :
                    plane[j,i] = np.mean(Weights[mask])
                    plane_var[j,i] = np.var(Weights[mask]) 
            rr_y += dR
        rr_x += dR
        
    if label == "mean" :
        return plane, Weights, dist_plane, thetas
    elif label == "var" :
        return plane_var, Weights, dist_plane, thetas
    else :
        return plane, plane_var, Weights, dist_plane, thetas




# Function that returns an array with elements in the unit disk.
# It is USELESS because the CreatePlane functions build arrays already inside the drive disk.
# KEEP for future documentation / development
# Use the MatrixMasked function for other purposes
def MatrixOnUnitDisk(Npixel, matrix) :
    
    # Initialize the Zernike2d class
    ZernikeM = ZF.Zernike2d(matrix)
    
    # Change of base of the matrix in the polar base
    matrix_polar = ZernikeM.r
    
    # I identify which distances are inside the unit disk (less than or equal to 1) through a mask
    matrix_mask = (matrix_polar <= 1)
    
    # New matrix with non-null values ​​only in the cells relative to the distances in the unit disk
    matrix_new = np.where(matrix_mask, matrix, 0)
    
    # Note: The last statement is equivalent to
    """
    matrix_new = np.zeros((Npixel,Npixel))
    for i in range(Npixel) :
        for j in range(Npixel) :
            if matrix_mask[i][j]:
                matrix_new[i][j] = matrix[i][j]
            else :
                matrix_new[i][j] = 0
    """
    
    return matrix_new




# Function that returns an array with non-null elements only in cells where the mask is True.
# It is needed in PlotMeanVariancePatch to plot a solid color for values below the chosen threshold.
def MatrixMasked(matrix, mask) :
    matrix_new = np.where(mask, matrix, 0)
    return matrix_new




# Function that returns the percentage of variances higher than a certain threshold of a patch on a drive disk.
# Also returns the relative mean and variance on the plane (required for PlotMeanVariancePatch input) 
def PercHigherVariance_Weights(label, Npixel, surf_a_obj, center, Dpp, threshold) :
    
    patch, mask = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
    rot_patch, _ = surf_a_obj.PatchReorientNew(patch, 1)
    
    ## To project the patch on the xy plane...
    z = surf_a_obj.FindOrigin(rot_patch)
    
    # Per ottenere il piano delle varianze (rispettivamente senza o con il piano delle medie)
    if label == "var" :
        plane_var, _, _, _ = CreatePlane_Weights("var", patch=rot_patch, z_c=z , Np=Npixel)
    else :
        plane, plane_var, _, _, _ = CreatePlane_Weights("", patch=rot_patch, z_c=z , Np=Npixel)
    
    # Inizializzo la classe Zernike2d
    ZernikeM = ZF.Zernike2d(plane_var)
    
    # Cambiamento di base della matrice nella base polare
    plane_var_polar = ZernikeM.r
    
    # Individuo quali distanze sono all'interno del disco unitario (cioè sono minori o uguali a 1) tramite una maschera
    polar_mask = (plane_var_polar <= 1)
    
    # Valori delle varianze nel disco unitario
    plane_var_masked = plane_var[polar_mask]
    
    # Numero di pixel all'interno del cerchio unitario
    Npixel_new = len(plane_var_masked)
    
    # Numero di varianze non nulle sopra la soglia (threshold) usando una maschera
    num_high_var = np.count_nonzero( plane_var_masked > threshold )
    
    perc = num_high_var / Npixel_new 
    
    if label == "var" :
        return plane_var, perc
    else :
        return plane, plane_var, perc




# PercHigherVariance version based on CreatePlane_Projections function. 
def PercHigherVariance_Projections(label, Npixel, surf_a_obj, center, Dpp, threshold) :
    
    patch, mask = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
    rot_patch, _ = surf_a_obj.PatchReorientNew(patch, 1)
    
    ## To project the patch on the xy plane...
    z = surf_a_obj.FindOrigin(rot_patch)
    
     # Per ottenere il piano delle varianze (rispettivamente senza o con il piano delle medie)
    if label == "var" :
        plane_var, _, _, _ = CreatePlane_Projections("var", patch=rot_patch, z_c=z , Np=Npixel)
    else :
        plane, plane_var, _, _, _ = CreatePlane_Projections("", patch=rot_patch, z_c=z , Np=Npixel)
    
    # Inizializzo la classe Zernike2d
    ZernikeM = ZF.Zernike2d(plane_var)
    
    # Cambiamento di base della matrice nella base polare
    plane_var_polar = ZernikeM.r
    
    # Individuo quali distanze sono all'interno del disco unitario (cioè sono minori o uguali a 1) tramite una maschera
    polar_mask = (plane_var_polar <= 1)
    
    # Valori delle varianze nel disco unitario
    plane_var_masked = plane_var[polar_mask]
    
    # Numero di pixel all'interno del cerchio unitario
    Npixel_new = len(plane_var_masked)
    
    # Numero di varianze non nulle sopra la soglia (threshold) usando una maschera
    num_high_var = np.count_nonzero( plane_var_masked > threshold )
    
    perc = num_high_var / Npixel_new 
    
    if label == "var" :
        return plane_var, perc
    else :
        return plane, plane_var, perc




# Function that produces the graph:
   # - The average of the distances that fall in each pixel as a function of the position (x, y)
   # - The variance of the distances that fall in each pixel as a function of the position (x, y)
#NOTE: both matrices have elements within the unit circle by construction

# References: 
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# https://matplotlib.org/stable/gallery/misc/rasterization_demo.html 
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

    fig, ax = mpl.subplots(nrows=1, ncols=2, figsize=(8,4), facecolor="white", dpi=200) # dpi=200 per compensare rasterized
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
            im = ax[col].pcolormesh(data, cmap=color_maps[col], rasterized=True) # senza rasterized il file è troppo grande
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
    



# Function that returns the index and value of the maximum and minimum of a portion of an array. 
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




# Function that returns:
# - the processed version of the average plane matrix, ie without empty areas.
# - the reconstruction according to Zernike of the latter.
# - the related Zernike complex cofficients. 
def PlaneRebuild(surf_a_obj, plane, ZOrder) :
    plane_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane) )
    zernike_env = ZF.Zernike2d(plane_proc)
    plane_recon, plane_coeff = zernike_env.ZernikeRecostruction(order=ZOrder, PLOT=0)
    return plane_proc, plane_recon, plane_coeff




# Function that produces the graph of:
# - the original version of the means matrix.
# - the processed version of the averages matrix.
# - the Zernike reconstructed version of the averages matrix.
def PlotPlaneRebuild(center, Dpp, Rs, ZOrder, pm, pm_p, pm_r, color_map, name) :
    
    if len(color_map) != 1 :
        color_map = "Greens"
   
    matrix = [ pm, pm_p, pm_r ]
    
    s0 = "Patch of center = {}, distance between points = {}, radius = {}, ZOrder = {}".format(center, Dpp, Rs, ZOrder)
    
    s1 = "Original Mean\n"
    s2 = "Processed Mean\n"
    s3 = "Reconstructed Mean\n"
        
    titles = [ s1, s2, s3 ]
        
    fig, ax = mpl.subplots(nrows=1, ncols=3, figsize=(12,4), dpi=200, facecolor="white")  # dpi=200 per compensare rasterized
    fig.suptitle(s0, fontsize="9")

    for row in range(1) :
        for col in range(3):
            data = matrix[col]
            ax[col].set_title(titles[col], fontsize="8")
            ax[col].set_xticks([])
            ax[col].set_yticks([])
            for side in ax[col].spines.keys():  # 'top', 'bottom', 'left', 'right'
                ax[col].spines[side].set_linewidth(0.30)
                ax[col].spines[side].set_color("black")
            im = ax[col].pcolormesh(data, cmap=color_map, rasterized=True)   # senza rasterized il file è troppo grande
            ax[col].axes.set_aspect("equal")  # necessario per correggere l'aspetto dei grafici senza colorbar
            ### colorbar ###
            """
            cb = mpl.colorbar(im, ax=ax[col], ticks=[])
            for side in cb.ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
                cb.ax.spines[side].set_linewidth(0.30)
                cb.ax.spines[side].set_color("black")
            """
            ### colorbar ###

    fig.tight_layout()
    
    if name != "" or name == "default" :
        if name == "default" :
            n = "Mean_Patch{}_Dpp{}_perc{}".format(center, Dpp, perc)
        else :
            n = name
        mpl.savefig("{}.pdf".format(n))
        print("The figure was generated.")

 

        
# Function that finds, in a group of points, the closest point to another point P.
# Will be used in GroupNearPoints to find the point of the contact zone closest to the center of mass of that zone.     
def PointNearPoint(points, P) :
    dist_points_P = np.sqrt( (points[:,0]-P[0])**2 + (points[:,1]-P[1])**2 + (points[:,2]-P[2])**2 )
    min_dist = np.amin(dist_points_P)
    point_near_P = int( np.where(dist_points_P == min_dist)[0] )
    return point_near_P




# Function that finds the contact zone between two proteins within a threshold distance (Daa).
# It, for each of the two proteins, returns:
# - the index of the point on the surface closest to the center of mass of that zone.
# - the matrix of the points in the contact area.
def GroupNearPoints(Daa, surf_a_obj, surf_b_obj) :
    
    prot_A = surf_a_obj.surface[:,:3]
    prot_B = surf_b_obj.surface[:,:3]
    print("Research of contact points. This step requires time...")
    patch_prot_a, patch_prot_b = SF.ContactPoints(prot_A, prot_B, Daa)
    print("Research complete.")
    
    cm_a = np.mean(patch_prot_a[:,:3], axis=0)
    print("CM of protein A group =", cm_a)
    center_a = PointNearPoint(patch_prot_a[:,:3], cm_a)
    print("Patch protein A: Center = {} with coord = {}".format(center_a, patch_prot_a[center_a,:3]))
    
    cm_b = np.mean(patch_prot_b[:,:3], axis=0)
    print("CM of protein B group =", cm_b)
    center_b = PointNearPoint(patch_prot_b[:,:3], cm_b)
    print("Patch protein B: Center = {} with coord = {}".format(center_b, patch_prot_b[center_b,:3]))
    
    return center_a, patch_prot_a[:,:3], center_b, patch_prot_b[:,:3]




# Function that generates the averages of two patches with opposite orientations (so as to be comparable)
# Each of the two patches belongs to a different protein.
# Its plans are generated with two different methods: CreatePlane_Weights and CreatePlane_Projections. 
def PatchesMethods(Npixel, surf_a_obj, c_a, surf_b_obj, c_b, Dpp) :
    
    patch_a, _ = surf_a_obj.BuildPatch(point_pos=c_a, Dmin=Dpp)
    rot_patch_a, rot_patch_nv_a = surf_a_obj.PatchReorientNew(patch_a, +1)
    z_pa = surf_a_obj.FindOrigin(rot_patch_a)
    plane_W_a, _, _, _ = CreatePlane_Weights("mean", patch=rot_patch_a, z_c=z_pa, Np=Npixel)
    plane_P_a, _, _, _ = CreatePlane_Projections("mean", patch=rot_patch_a, z_c=z_pa, Np=Npixel)
    
    patch_b, _ = surf_b_obj.BuildPatch(point_pos=c_b, Dmin=Dpp)
    rot_patch_b, rot_patch_nv_b = surf_b_obj.PatchReorientNew(patch_b, -1)
    z_pb = surf_b_obj.FindOrigin(rot_patch_b)
    plane_W_b, _, _, _ = CreatePlane_Weights("mean", patch=rot_patch_b, z_c=z_pb, Np=Npixel)
    plane_P_b, _, _, _ = CreatePlane_Projections("mean", patch=rot_patch_b, z_c=z_pb, Np=Npixel)
    
    return plane_W_a, plane_P_a, plane_W_b, plane_P_b



# Function that calculates, more quickly than PlaneRebuild, the Zernike coefficients and their modules
# To have the right coefficients the plan must be the one processed NOT the original one
# The processed plane is the one that is flooded and without empty pixels 
def ZernikeCoeff(ZOrder, surf_a_obj, plane) :
    plane_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane) )
    zernike_env = ZF.Zernike2d(plane_proc)
    coeff = zernike_env.ZernikeDecomposition(order=ZOrder)  # coeff is a list
    coeff_inv = np.absolute(coeff)
    return plane_proc, coeff, coeff_inv




# Function that calculates the difference between two lists of Zernikei coefficient modules 
def ZernikeCoeff_Distance(ZOrder, surf_a_obj, plane_1, surf_b_obj, plane_2) :
    _, _, c_inv_1 = ZernikeCoeff(ZOrder, surf_a_obj, plane_1)
    _, _, c_inv_2 = ZernikeCoeff(ZOrder, surf_b_obj, plane_2)
    c_inv_diff = np.sqrt( sum( (c_inv_1[:]-c_inv_2[:])**2 ) )
    return c_inv_diff




# Function that compares the average of a patch produced with two different methods: CreatePlane_Weights and CreatePlane_Projections. 
def PlotPatchesComparison(obj_name, Npixel, Rs, p_W, p_P, center, Dpp, Daa, color_map, name) :
    
    if obj_name == "" :
        obj_name = "Unknown"

    if len(color_map) != 1 :
        color_map = "Greens"
   
    matrix = [ p_W, p_P ]
    
    s0 = "Protein {} with threshold distance = {}\nPatch of center = {}, distance between points = {}, radius = {}".format(obj_name, Daa, center, Dpp, Rs)
    
    s1 = "Mean with Weights Method in function of position\n\n"
    s2 = "Mean with Projections Method in function of position\n\n"
    
    titles = [ s1, s2 ]
        
    fig, ax = mpl.subplots(nrows=1, ncols=2, figsize=(8,4), dpi=200, facecolor="white")  # dpi=200 per compensare rasterized
    fig.suptitle(s0, fontsize="9")

    for row in range(1) :
        for col in range(2):
            data = matrix[col]
            ax[col].set_title(titles[col], fontsize="8")
            ax[col].set_xlabel("x", fontsize="8")
            ax[col].set_ylabel("y", fontsize="8")
            ax[col].tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in ax[col].spines.keys():  # 'top', 'bottom', 'left', 'right'
                ax[col].spines[side].set_linewidth(0.30)
                ax[col].spines[side].set_color("black")
            im = ax[col].pcolormesh(data, cmap=color_map, rasterized=True)   # senza rasterized il file è troppo grande
            #ax[col].axes.set_aspect("equal")  # necessario per correggere l'aspetto dei grafici senza colorbar
            ### colorbar ###
            ticks_list = [np.amin(data), np.amax(data)]
            cb = mpl.colorbar(im, ax=ax[col], ticks=ticks_list)
            cb.ax.tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in cb.ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
                cb.ax.spines[side].set_linewidth(0.30)
                cb.ax.spines[side].set_color("black")
            ### colorbar ###

    fig.tight_layout()
    
    if name != "" or name == "default" :
        if name == "default" :
            n = "Mean_Patch{}.pdf".format(center, Dpp)
        else :
            n = name
        mpl.savefig("{}".format(n))
        print("The figure was generated.")
