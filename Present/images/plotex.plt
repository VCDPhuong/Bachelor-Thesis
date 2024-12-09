set term pdfcairo enhanced color
set output  "Experiment.pdf" 
set xlabel "Energy (eV)"
set ylabel "Absorption Spectrum"
plot "plot-data.csv" w l lw 2 notitle, "plot-data.csv" w labels center offset -0.2,0.7 notitle 