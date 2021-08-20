**Autore**: Alessandro Giudice    
**Collaboratore**: Samuel Santhosh Gomez  

# Percentuale di non funzionalità di una patch in disco unitario
Di seguito riporto la procedura per calcolare la percentuale di non funzionalità di una specifica zona di una superficie proteica in 3D.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, che è visibile in appendice oppure <a href="https://github.com/agds93/percentage_non_functionality/tree/main/code" target="_blank">qui</a>.   

## Selezione di una patch
L'intera superficie proteica studiata è visibile in Figura 0.

<p align="center"><img src="img/entire_protein.png" width=600px></p>
<p align="center"><i>Figura 0</i>: L'intera superficie proteica in 3D.</p>

Una patch, come quella in Figura 1, è un gruppo di punti di una superficie 3D. Tali punti sono selezionati come punti della patch se essi hanno una distanza reciproca non superiore al valore soglia `Dpp`, e se sono contenuti in una sfera avente come centro l'indice di un punto della superficie `center` e un raggio `Rs`.

<p align="center"><img src="img/Patch_Point5000.png" width=440px><img src="img/Patch_Point19841.png" width=440px></p>
<p align="center"><i>Figura 1</i>: Patch del punto 5000 (sinistra) e patch del punto 19841 (destra) della superficie.</p>

La patch selezionata deve essere ruotata per essere perpendicolare al piano xy, poi è inglobata in un cono come in Figura 2. Tale cono è posto lungo l'asse z, con origine nel punto C=(0,0,`z`), in modo che l'angolo massimo tra l'asse perpendicolare e la secante che connette C a un punto della superficie (o della patch) sia uguale a `theta_max = 45`.

<p align="center"><img src="img/Cone_Point5000.png" width=600px></p>
<p align="center"><i>Figura 2</i>: La patch (rosso) del punto 5000 all'interno del cono (blu).</p>

## Creazione del piano di fit
Ogni punto della patch viene proiettato su una griglia quadrata 2D di lato `Npixel`, in cui ogni cella è un pixel. All'interno di ogni pixel è presente il valore della media o della varianza delle distanze tra i relativi punti della patch e l'origine C del cono. Di conseguenza le possibili griglie da creare sono due:

* la matrice delle medie delle distanze in ogni pixel, come nella parte sinistra delle Figure 3-4.

* la matrice delle varianze delle distanze in ogni pixel, come nella parte destra delle Figure 3-4.

Le distanze utilizzate per ogni matrice nelle Figure 3-4 sono solo quelle contenute in un disco unitario (distanze minori o uguali a uno).  
Per creare il piano della media e della varianza di tali distanze ci sono due metodi.  
Il primo metodo (funzione `CreatePlane_Weights`) costruisce una griglia in cui il valore (media o varianza) di ogni pixel si basa sulle distanze tra i punti della patch e il punto C. La distanza relativa ad un punto della patch finisce in un pixel se il punto si trova ortogonalmente sopra tale pixel. Due esempi di tale metodo sono visibili in Figura 3. Le parti alta e bassa della figura sono riferite rispettivamente alla patch con `center = 5000` e alla patch con `center = 19841`. Questo metodo può essere chiamato metodo *Weights*.

<p align="center">
<img src="img/Point_5000_Weights.png" width=700px>
<img src="img/Point_19841_Weights.png" width=700px>
</p>
<p align="center"><i>Figura 3</i>: Media e varianza di due patch (una per riga) prodotte con il primo metodo (metodo <i>Weights</i>).</p>

Nel secondo metodo (funzione `CreatePlane_Projections`) la griglia viene costruita in modo che ogni pixel abbia un valore (media o varianza) basato sulle distanze tra i punti della patch e il punto C. A differenza del primo metodo, la distanza relativa ad un punto della patch finisce in un pixel se il segmento che congiunge un punto della patch e il punto C intercetta tale pixel. Gli stessi esempi di Figura 3 prodotti con tale metodo sono visibili in Figura 4. Questo metodo può essere chiamato metodo *Projections*.

<p align="center">
<img src="img/Point_5000_Projections.png" width=700px>
<img src="img/Point_19841_Projections.png" width=700px>
</p>
<p align="center"><i>Figura 4</i>: Media e varianza di due patch (una per riga) prodotte con il secondo metodo (metodo <i>Projections</i>).</p>

## Percentuale di non funzionalità
La percentuale di non funzionalità `perc` di una patch coincide con la percentuale di pixels della matrice che contengono una varianza superiore ad una soglia `threshold`. Il valore trovato di `perc` e il valore scelto per `threshold` è riportato nel titolo della parte destra dei grafici delle Figure 3-4. Inoltre il valore della soglia è indicato anche sulla relativa barra colorata di tali figure. Per ogni pixel, se la varianza è inferiore a tale soglia viene mostrato un colore uniforme (patch con `center = 5000` in Figura 3-4), in caso contrario viene visualizzato un colore più o meno scuro per un valore alto o basso della varianza (patch con `center = 19841` in Figura 3-4).  
I valori di `perc` sono calcolati con le funzioni `PercHigherVariance_Weights` e `PercHigherVariance_Projections`. Tali valori per ogni punto della superficie sono visibili in Figura 4 e Figura 5 rispettivamente per primo e secondo metodo.

<p align="center"><img src="img/all_perc.png" width=800px></p>
<p align="center"><i>Figura 4</i>: Percentuale di non funzionalità con il metodo <i>Weights</i> per ogni punto della superficie.</p>
<p align="center"><img src="img/all_perc_projections.png" width=800px></p>
<p align="center"><i>Figura 5</i>: Percentuale di non funzionalità con il metodo <i>Projections</i> per ogni punto della superficie.</p>

Come mostrato in Figura 6-7, il secondo metodo produce piani di fit con una percentuale di non-funzionalità generalmente più bassa rispetto al primo metodo.

<p align="center"><img src="img/hist_01.png" width=800px></p>
<p align="center"><i>Figura 6</i>: Istogramma della percentuale di non funzionalità con il metodo <i>Weights</i>.</p>
<p align="center"><img src="img/hist_02.png" width=800px></p>
<p align="center"><i>Figura 7</i>: Istogramma della percentuale di non funzionalità con il metodo <i>Projections</i>.</p>

## Appendice
### Librerie e moduli
Il codice scritto è stato eseguito con <a href="https://jupyterlab.readthedocs.io/en/stable/" target="_blank">JupyterLab</a> utilizzando `python 3.8`.  
I moduli python usati, installati tramite 
<a href="https://pip.pypa.io/en/stable/" target="_blank">pip</a> (compreso `jupyterlab`), sono elencati sotto.
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
Il modulo `mayavi`, in particolare `mlab`, è necessario per visualizzare le superfici 3D in una finestra Qt, così da produrre la Figura 0-1-2.  
Mentre le librerie di base sono
```python
sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF
```
scritte da <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.

### Parametri
I valori dei parametri usati per selezionare una patch e creare il piano di fit sono
```python
Npixel = 25    # il lato del piano in pixel
Dpp = 0.5      # la distanza tra i punti della stessa patch
Rs = 6         # il raggio della sfera che include la patch
threshold = 5  # valore soglia per stabilire se la varianza è alta
```
I valori sono in ångström, tranne `Npixel`.  
L'indice del punto centrale della patch usato per Figura 1, Figura 2, e per i grafici nella parte alta della Figura 3-4 è
```python
center = 5000
```
invece quello usato per i grafici nella parte bassa della Figura 3-4 è
```python
center = 19841
```
### Caricare la superficie proteica
Una volta scelta la proteina da studiare bisogna scaricare il relativo file *.pdb* dalla
<a href="https://www.rcsb.org/" target="_blank">Protein Data Bank</a>
da cui creare un file *.dms* (per esempio con il tool
<a href="https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html" target="_blank">dms</a>
), contenente una serie di dati su atomi e punti della superficie.
```python
surf_name_a = "./data/4bs2_RRM2.dms"
surf_a_ = pd.read_csv(surf_name_a)  
l_a = len(surf_a_["x"])
print("Npoints", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
```
dove `surf_name_a` è il percorso del file *.dms* usato, disponibile <a href="data/4bs2_RRM2.dms" target="_blank">qui</a>.  
La matrice relativa all'intera superficie `surf_a` deve essere inizializzato come oggetto della classe `Surface`: 
```python
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
Dopo aver caricato i punti della superficie completa, per graficare l'intera superficie proteica come in Figura 0 si usa
```python
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
### Selezione di una patch
Una patch è costruita in base ai parametri scelti tramite
```python
patch, _ = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
```
Per produrre il grafico in Figura 1, con `center = 5000`, si usa:
```python
res1, c = SF.ConcatenateFigPlots([patch[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
Per essere utilizzabile, la patch in Figura 1 deve essere ruotata (`patch` diventa `rot_patch`) in modo che sia perpendicolare al piano xy. Questo si può fare con
```python 
rot_patch, rot_patch_nv = surf_a_obj.PatchReorientNew(patch, +1)
```
dove il parametro `+1` (`-1`) indica che i versori normali `rot_patch_nv` sono rivolti verso l'alto (basso).  
Per trovare la quota `z` dell'origine C del cono che ingloba la patch si usa
```python
z = surf_a_obj.FindOrigin(rot_patch, 0)
```
dove se si sostituisce `0` con `1` si produce il grafico della patch inglobata dentro il cono. La Figura 2 si ottiene con `center = 5000`.
### Creazione del piano di fit
Per creare la matrice della media `plane` e delle varianza `plane_var` con il primo metodo si utilizza la funzione seguente. L'input è formato da:
* una label che determina se restituire la matrice di media e/o la varianza:
    * "mean" restituisce solo i valori medi
    * "var" produce solo le varianze
    * "" o altro fornisce valori medi e varianze
* la patch ruotata.
* la quota `z` del punto C.
* il numero di pixel per lato della griglia `Npixel`.
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
Invece per creare la matrice della media `plane` e delle varianza `plane_var` con il secondo metodo si utilizza la funzione seguente. Gli input sono gli stessi della funzione utilizzata per il primo metodo.
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
### Percentuale di non funzionalità
Per calcolare la percentuale `perc` di varianze più alte di una certa soglia con il primo metodo si usa la funzione seguente. L'input è formato da:
* una label che determina se restituire, oltre a `perc`, la matrice di media e/o varianza:
    * "mean" restituisce solo i valori medi
    * "var" produce solo le varianze
    * "" o altro fornisce valori medi e varianze
* il numero di pixel per lato della griglia `Npixel`.
* l'oggetto superficie `surf_a_obj`.
* l'indice `center` del punto della superficie scelto come centro della patch.
* la distanza `Dpp` tra i punti della patch.
* la soglia `threshold` scelta.

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
Invece per calcolare la percentuale `perc` di varianze più alte di una certa soglia con il secondo metodo si usa la funzione seguente. Gli input sono gli stessi della funzione utilizzata per il primo metodo, infatti l'unica differenza è l'utilizzo di `CreatePlane_Projections` invece di `CreatePlane_Weights`.
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
Per calcolare la `perc` di ogni punto della superficie con il primo metodo si usa
```python
limit = l_a
step = 1
points_list = np.arange(0,limit,step)
perc = np.zeros((len(points_list)))
for i in range(len(points_list)) :
    _, _, perc[i] = PercHigherVariance_Weights(Npixel, Rs, surf_a_obj, points_list[i], Dpp, threshold)  
print("Number of patches =",len(perc))
```
Per calcolare la `perc` di ogni punto della superficie con il secondo metodo si usa
```python
limit = l_a
step = 1
points_list = np.arange(0,limit,step)
perc = np.zeros((len(points_list)))
for i in range(len(points_list)) :
    _, _, perc[i] = PercHigherVariance_Projections(Npixel, Rs, surf_a_obj, points_list[i], Dpp, threshold)
print("Number of patches =",len(perc))
```
Tale codice impiega molto tempo per essere eseguito quindi i risultati sono salvati su un file:
* i valori trovati con `PercHigherVariance_Weights` sono visibili <a href="/data/all_perc.txt" target="_blank">qui</a>.
* i valori trovati con `PercHigherVariance_Projections` sono visibili <a href="/data/all_perc_projections.txt" target="_blank">qui</a>.
```python
with open("all_perc.txt", "w") as file0 :
    for i in range(len(points_list)) :
        file0.write("{}\t{}\n".format(points_list[i],perc[i]))
```
così è possibile caricare indici delle patch e relative percentuali direttamente dal file.
```python
points_list = np.loadtxt("./risultati/all_perc.txt", usecols=0, unpack=True)
perc = np.loadtxt("./risultati/all_perc.txt", usecols=1, unpack=True)
```
Per trovare i massimi e i minimi di `perc` in uno spefico intervallo si usa 
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
### Utilità
Segue la funzione che restituisce una matrice con elementi non nulli solo nelle celle in cui la maschera ha valore True. Essa è necessaria in `PlotMeanVariancePatch` per graficare un colore uniforme per i valori della varianza sotto la soglia scelta.
```python
def MatrixMasked(matrix, mask) :
    matrix_new = np.where(mask, matrix, 0)
    return matrix_new
```
### Altri Grafici
La seguente funzione fornisce i grafici delle Figure 3-4. Gli input sono:
* l'indice del punto centrale della patch `center`.
* la distanza tra i punti della patch `Dpp`.
* il raggio `Rs` della sfera che include la patch.
* il valore trovato di `perc`.
* il valore di soglia scelto.
* la matrice della media.
* la matrice della varianza.
* le mappe dei <a href="https://matplotlib.org/stable/tutorials/colors/colormaps.html" target="_blank">colori</a> da utilizzare.
* il nome del file di output con opportuna estensione.
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
Il grafico in Figura 4-5 viene prodotto da
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
Il grafico in Figura 6-7 viene prodotto da
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