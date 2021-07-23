Di seguito riporto la procedura per calcolare la percentuale di non funzionalità di una specifica zona di una superficie proteica in 3D.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, visibile in appendice.   

## Selezione di una patch
Una patch, come quella in Figura 1, è un gruppo di punti di una superficie 3D.
<figure class="image">
  <img src="img/Patch_Point5000.png" width=600px class="center">
  <center><figcaption><i>Figura 1</i>: Una possibile patch della superficie.</figcaption></center>
</figure>
Tali punti sono selezionati come punti della patch se essi hanno una distanza reciproca non superiore al valore soglia `Dpp`, e se sono contenuti in una sfera avente come centro l'indice di un punto dell'intera superficie `center` e un raggio `Rs`.  
Poi la patch selezionata deve essere inglobata in un cono come in Figura 2.
<figure class="image">
  <img src="img/Cone_Point5000.png" width=600px class="center">
  <center><figcaption><i>Figura 2</i>: Patch (rosso) all'interno del cono (blu).</figcaption></center>
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
  <img src="img/Point_19841_Projections.png" width=600px class="center">
  <center><figcaption><i>Figura 3</i>: Media e varianza di due diverse patch prodotte con il primo metodo.</figcaption></center>
</figure>
<figure class="image">
  <img src="img/Point_5000_Weigths.png" width=600px class="center">
  <img src="img/Point_19841_Projections.png" width=600px class="center">
  <center><figcaption><i>Figura 4</i>: Media e varianza di due diverse patch prodotte con il secondo metodo.</figcaption></center>
</figure>
