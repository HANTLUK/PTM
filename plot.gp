set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"

set decimalsign "." # für den input
#               Style
# !!!___________________________ !!!
set grid xtics
set grid ytics
# set grid mxtics
# set grid mytics
set style line 80 linetype 1 linecolor rgb "#888888"
set style line 81 linetype 1 linecolor rgb "#808080" linewidth 0.5
set border back linestyle 80
set grid back linestyle 81
set xtics textcolor rgb "#808080"
set ytics textcolor rgb "#808080"
set y2tics textcolor rgb "#808080"

LINECOLORS = 'red pink'
LINEWIDTHS = '1.0  4.0   0.0   0.0     0.0'
myLinecolor(i) = word(LINECOLORS, i)
myLinewidth(i) = real(LINEWIDTHS, i)

set xlabel "Number of Qubits $n$"
set ylabel "Execution Time $s$"

set xrange [0.5:7.5]
set yrange [0.0000001:1000]

set logscale y
set format y "$10^{%T}$"

set nokey

# !!!_______________________________________________!!!
#		C-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/C-PTM.tex"

set title "C-PTM"
plot  'TabsPTM/C-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 2 ps 0.2 lc "#1188ff" ,\
'TabsPTM/C-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 7 lw 2 ps 0.2 lc "#1188ff"
# !!!_______________________________________________!!!
#		L-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/L-PTM.tex"

set title "L-PTM"
plot  'TabsPTM/L-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 2 ps 0.2 lc "#1188ff" ,\
'TabsPTM/L-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 7 lw 2 ps 0.2 lc "#1188ff"
# !!!_______________________________________________!!!
#		M-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/M-PTM.tex"

set title "M-PTM"
plot  'TabsPTM/M-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 2 ps 0.2 lc "#1188ff" ,\
'TabsPTM/M-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 7 lw 2 ps 0.2 lc "#1188ff"
# !!!_______________________________________________!!!
#		AC-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/AC-PTM.tex"

set title "AC-PTM"
plot  'TabsPTM/AC-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 2 ps 0.2 lc "#1188ff" ,\
'TabsPTM/AC-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 7 lw 2 ps 0.2 lc "#1188ff"

set xrange [0.5:9.5]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Kraus-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/Kraus-PTM.tex"

set title "Kraus-PTM"

plot 'TabsPTM/Kraus-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" ps 0.2 ,\
'TabsPTM/Kraus-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" lw 2 ,\
'TabsPTM/Kraus-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Kraus-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" dashtype 4 lw 2 ,\
'TabsPTM/Kraus-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" ps 0.2 ,\
'TabsPTM/Kraus-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" lw 2 ,\
'TabsPTM/Kraus-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Kraus-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" dashtype 4 lw 2


set xrange [0.5:9.5]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Choi-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/Choi-PTM.tex"

set title "Choi-PTM"

plot 'TabsPTM/Choi-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" ps 0.2 ,\
'TabsPTM/Choi-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" lw 2 ,\
'TabsPTM/Choi-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Choi-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" dashtype 4 lw 2 ,\
'TabsPTM/Choi-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" ps 0.2 ,\
'TabsPTM/Choi-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" lw 2 ,\
'TabsPTM/Choi-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Choi-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" dashtype 4 lw 2 ,\
'TabsPTM/Choi-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" ps 0.2 ,\
'TabsPTM/Choi-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" lw 2 ,\
'TabsPTM/Choi-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Choi-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" dashtype 4 lw 2

set xrange [0.5:9.5]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Chi-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/Chi-PTM.tex"

set title "Chi-PTM"

plot 'TabsPTM/Chi-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" ps 0.2 ,\
'TabsPTM/Chi-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" lw 2 ,\
'TabsPTM/Chi-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Chi-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" dashtype 4 lw 2 ,\
'TabsPTM/Chi-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" ps 0.2 ,\
'TabsPTM/Chi-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" lw 2 ,\
'TabsPTM/Chi-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Chi-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" dashtype 4 lw 2 ,\
'TabsPTM/Chi-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" ps 0.2 ,\
'TabsPTM/Chi-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" lw 2 ,\
'TabsPTM/Chi-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#33ff44" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Chi-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#11bb22" dashtype 4 lw 2


set xrange [0.5:9.5]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Can-PTM

set term cairolatex pdf size 7cm,4.5cm color colortext font ",8"
set output "FigsPTM/Can-PTM.tex"

set title "Can-PTM"

plot 'TabsPTM/Can-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" ps 0.2 ,\
'TabsPTM/Can-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" lw 2 ,\
'TabsPTM/Can-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Can-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#1144bb" dashtype 4 lw 2 ,\
'TabsPTM/Can-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" ps 0.2 ,\
'TabsPTM/Can-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" lw 2 ,\
'TabsPTM/Can-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff3344" dashtype 4 lw 2 ps 0.2 ,\
'TabsPTM/Can-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#bb1122" dashtype 4 lw 2
