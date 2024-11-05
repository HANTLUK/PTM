set term cairolatex pdf size 8cm, 8cm color colortext font ",8"

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

# "#0072B2, #009E73, #D55E00, #CC79A7, #F0E442, #56B4E9"

set xlabel "Number of Qubits $n$"
set ylabel "Execution Time $s$"

set xrange [0.9:7.1]
set yrange [0.0000001:1000]

set logscale y
set format y "$10^{%T}$"

set nokey

# !!!_______________________________________________!!!
#		C-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/C-PTM.tex"

set title "\\texttt\{C-PTM\}"
plot  'TabsPTM/C-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 4 ps 0.6 lc "#1188ff" ,\
'TabsPTM/C-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 4 lw 4 ps 0.6 lc "#1188ff"
# !!!_______________________________________________!!!
#		L-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/L-PTM.tex"

set title "\\texttt\{L-PTM\}"
plot  'TabsPTM/L-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 4 ps 0.6 lc "#1188ff" ,\
'TabsPTM/L-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 4 lw 4 ps 0.6 lc "#1188ff"
# !!!_______________________________________________!!!
#		M-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/M-PTM.tex"

set title "\\texttt\{M-PTM\}"
plot  'TabsPTM/M-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 4 ps 0.6 lc "#1188ff" ,\
'TabsPTM/M-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 4 lw 4 ps 0.6 lc "#1188ff"
# !!!_______________________________________________!!!
#		AC-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/AC-PTM.tex"

set title "\\texttt\{AC-PTM\}"
plot  'TabsPTM/AC-PTM_randDense_TPD.dat' using 1:2 with linespoints lw 4 ps 0.6 lc "#1188ff" ,\
'TabsPTM/AC-PTM_randDiag_TPD.dat' using 1:2 with linespoints dashtype 4 lw 4 ps 0.6 lc "#1188ff"

set xrange [0.9:9.1]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Kraus-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/Kraus-PTM.tex"

set title "\\texttt\{Kraus-PTM\}"

plot 'TabsPTM/Kraus-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" lw 4 ps 0.6 ,\
'TabsPTM/Kraus-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" lw 4 ,\
'TabsPTM/Kraus-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Kraus-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" dashtype 4 lw 4 ,\
'TabsPTM/Kraus-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" lw 4 ps 0.6 ,\
'TabsPTM/Kraus-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" lw 4 ,\
'TabsPTM/Kraus-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Kraus-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" dashtype 4 lw 4


set xrange [0.9:9.1]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Choi-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/Choi-PTM.tex"

set title "\\texttt\{Choi-PTM\}"

plot 'TabsPTM/Choi-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" lw 4 ,\
'TabsPTM/Choi-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" dashtype 4 lw 4 ,\
'TabsPTM/Choi-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" lw 4 ,\
'TabsPTM/Choi-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" dashtype 4 lw 4 ,\
'TabsPTM/Choi-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" lw 4 ,\
'TabsPTM/Choi-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Choi-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" dashtype 4 lw 4

set xrange [0.9:9.1]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Chi-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/Chi-PTM.tex"

set title "\\texttt\{Chi-PTM\}"

plot 'TabsPTM/Chi-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" lw 4 ,\
'TabsPTM/Chi-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" dashtype 4 lw 4 ,\
'TabsPTM/Chi-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" lw 4 ,\
'TabsPTM/Chi-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" dashtype 4 lw 4 ,\
'TabsPTM/Chi-PTM_randDense_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDense_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" lw 4 ,\
'TabsPTM/Chi-PTM_randDiag_HY.dat' using 1:2 with linespoints lc rgb "#11ff88" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Chi-PTM_randDiag_HY_ex.dat' using 1:2 with lines lc rgb "#02aa35" dashtype 4 lw 4


set xrange [0.9:9.1]
set yrange [0.00001:1000000]

# !!!_______________________________________________!!!
#		Can-PTM

set term cairolatex pdf size 8cm, 8cm color colortext font ",8"
set output "FigsPTM/Can-PTM.tex"

set title "\\texttt\{Can-PTM\}"

plot 'TabsPTM/Can-PTM_randDense_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" lw 4 ps 0.6 ,\
'TabsPTM/Can-PTM_randDense_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" lw 4 ,\
'TabsPTM/Can-PTM_randDiag_QK.dat' using 1:2 with linespoints lc rgb "#1188ff" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Can-PTM_randDiag_QK_ex.dat' using 1:3 with lines lc rgb "#0235aa" dashtype 4 lw 4 ,\
'TabsPTM/Can-PTM_randDense_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" lw 4 ps 0.6 ,\
'TabsPTM/Can-PTM_randDense_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" lw 4 ,\
'TabsPTM/Can-PTM_randDiag_TPD.dat' using 1:2 with linespoints lc rgb "#ff1188" dashtype 4 lw 4 ps 0.6 ,\
'TabsPTM/Can-PTM_randDiag_TPD_ex.dat' using 1:2 with lines lc rgb "#aa0235" dashtype 4 lw 4
