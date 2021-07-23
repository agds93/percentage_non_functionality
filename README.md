Di seguito riporto la procedura per calcolare la percentuale di non funzionalità di una specifica zona di una superficie proteica in 3D.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, visibile in appendice.   

## Selezione di una patch
Una patch, come quella in Figura 1, è un gruppo di punti di una superficie 3D.
<figure class="image">
  <img src="img/Patch_Point5000.png" width=600px class="center">
  <figcaption><i>Figura 1</i>: Una possibile patch della superficie.</figcaption>
</figure>

Tali punti sono selezionati come punti della patch se essi hanno una distanza reciproca non superiore al valore soglia `Dpp`, e se sono contenuti in una sfera avente come centro l'indice di un punto dell'intera superficie `center` e un raggio `Rs`.  
Poi la patch selezionata deve essere inglobata in un cono come in Figura 2.
<figure class="image">
  <img src="img/Cone_Point5000.png" width=600px class="center">
  <center><figcaption><i>Figura 2</i>: Patch (rosso) all'interno del cono (blu).  </figcaption></center>
</figure>

Tale cono è posto lungo l'asse z, con origine nel punto $C$=$(0,0,$ `z`$)$, in modo che l'angolo massimo tra l'asse perpendicolare e la secante che connette $C$ a un punto della superficie (o della patch) sia uguale a `theta_max = 45`.  


Come anticipato sopra, un piano di fit per sarà rappresentato da una matrice quadrata 2D di lato `Npixel`. Le matrici, cioè i piani di fit, producibili sono due:

* la matrice in cui in ogni elemento è la media delle distanze contenute nel relativo pixel.

* la matrice in cui in ogni elemento è la varianza delle distanze contenute nel relativo pixel.

## Creazione del piano di fit
Ogni punto della patch viene proiettato su una griglia quadrata 2D di lato `Npixel`, in cui ogni cella è un pixel. All'interno di ogni pixel ci possono essere il valore della media o della varianza delle distanze tra i relativi punti della patch e l'origine $C$ del cono. Di conseguenza i possibili piani da creare sono due:

* la matrice delle medie delle distanze in ogni pixel, come nella parte sinistra delle Figure 3-4.

* la matrice delle varianze delle distanze in ogni pixel, come nella parte destra delle Figure 3-4.

Le distanze utilizzate per ogni matrice nelle Figure 3-4 sono solo quelle contenute in un disco unitario (distanze minori o uguali a uno).

<figure class="image">
  <img src="img/Point_5000_Weigths.png" width=600px class="center">
  <img src="img/Point_19841_Weigths.png" width=600px class="center">
  <center><figcaption><i>Figura 3</i>: Media e varianza di due diverse patch prodotte con il primo metodo.</figcaption></center>
</figure>

<figure class="image">
  <img src="img/Point_5000_Projections.png" width=600px class="center">
  <img src="img/Point_19841_Projections.png" width=600px class="center">
  <center><figcaption><i>Figura 4</i>: Media e varianza di due diverse patch prodotte con il secondo metodo.</figcaption></center>
</figure>

## Appendice
### Librerie
Il codice scritto è stato eseguito con <a href="https://jupyterlab.readthedocs.io/en/stable/" target="_blank">JupyterLab</a> utilizzando almeno `python 3.8.10`.  
I moduli python usati, compreso `jupyterlab`, che sono stati installati tramite 
<a href="https://pip.pypa.io/en/stable/" target="_blank">pip</a>, sono elencati sotto.
```{.python}
import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import pandas as pd
```
```{.python}
from mayavi import mlab
```
in particolare `mlab` è necessario per visualizzare le superfici 3D tramite una finestra Qt, così da produrre le Figure 1-2.     
Mentre le librerie di base sono
```{.python}
sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF
```
scritte da <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.

### Parametri
I valori dei parametri usati per selezionare e fittare una patch sono
```{.python}
Npixel = 25    # il lato del piano in pixel
Dpp = 0.5      # la distanza tra i punti della stessa patch
Rs = 6         # il raggio della sfera che include la patch
threshold = 5  # valore soglia per stabilire se la varianza è alta (in ångström)
```
Gli indici della superficie usati per i grafici mostrati sono
```{.python}
center = 5000
```
```{.python}
center = 19841
```
### Caricare una superficie proteica
Una volta scelta la proteina da studiare bisogna scaricare il relativo file *.pdb* dalla
<a href="https://www.rcsb.org/" target="_blank">Protein Data Bank</a>
da cui creare un file *.dms* (per esempio con il tool
<a href="https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html" target="_blank">dms</a>
), contenente una serie di dati su atomi e punti della superficie.
```{.python}
surf_name_a = "./data/4bs2_RRM2.dms"
surf_a_ = pd.read_csv(surf_name_a)  
l_a = len(surf_a_["x"])
print("Npoints", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
```
dove `surf_name_a` è il percorso del file *.dms*.
La matrice relativa all'intera superficie `surf_a` deve essere inizializzato come oggetto della classe `Surface`: 
```{.python}
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
### Selezione di una patch
Dopo aver caricato la superficie completa, una patch è costruita in base ai parametri scelti tramite
```{.python}
patch, _ = surf_a_obj.BuildPatch(point_pos=center, Dmin=Dpp)
```
Per produrre il grafico in Figura 1, con `center = 5000`, si usa:
```{.python}
res1, c = SF.ConcatenateFigPlots([patch[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
Per essere utilizzabile, la patch deve essere ruotata (`patch` diventa `rot_patch`) in modo che sia perpendicolare al piano $xy$. Questo si può fare con
```{.python}  
rot_patch, rot_patch_nv = surf_a_obj.PatchReorientNew(patch, +1)
```
dove il parametro `+1` indica che i versori normali `rot_patch_nv` sono rivolti verso l'alto.  
Per trovare la quota `z` dell'origine $C$ del cono che ingloba la patch si usa
```{.python}
z = surf_a_obj.FindOrigin(rot_patch, 0)
```
dove se si sostituisce `0` con `1` si produce la Figura 2 con `center = 5000`.
