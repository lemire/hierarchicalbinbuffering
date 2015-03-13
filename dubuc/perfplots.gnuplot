set terminal postscript eps mono 25
set nokey

set ylabel "IA-32 floating pt. instr."
set output "fpins-varying-range128-sum.eps"
set yrange [250000:350000]
set xlabel "range length"
plot "rangesums_N=2_b=128_fpins.dat" with points
set output
!epstopdf  fpins-varying-range128-sum.eps
!mv fpins-varying-range128-sum.pdf ../../LaTeX/Ola/LaTeX

set ylabel "IA-32 floating pt. instr."
set output "fpins-varying-range128-fm.eps"
set yrange [250000:350000]
set xlabel "range length"
plot  "firstmoments_N=2_b=128_fpins.dat" with points
set output
!epstopdf  fpins-varying-range128-fm.eps
!mv fpins-varying-range128-fm.pdf ../../LaTeX/Ola/LaTeX


set ylabel "time (seconds)"
set output "construction-b.eps"
set logscale x
set xlabel "b"
set yrange [1450:1550]
plot "construction_N=2.dat" with points pointtype 13 pointsize 2
set output
!epstopdf  construction-b.eps
!mv construction-b.pdf ../../LaTeX/Ola/LaTeX
set yrange [*:*]


set output "time-varying-b.eps"
set logscale xy
set xrange[10:*]
plot "query_N=2.dat" with points pointtype 13 pointsize 2, x/(5000 * log(x))
set output
!epstopdf  time-varying-b.eps
!mv time-varying-b.pdf ../../LaTeX/Ola/LaTeX
set xrang [*:*]
set nologscale xy

set output "construction-N-smallb.eps"
set xlabel "N"
plot "construction_b=128.dat" with points pointtype 13 pointsize 2, 268*x+450
set output
!epstopdf  construction-N-smallb.eps
!mv construction-N-smallb.pdf ../../LaTeX/Ola/LaTeX

set output "construction-N-bigb.eps"
plot "construction_b=32k.dat" with points pointtype 13 pointsize 2,  268*x+450
set output
!epstopdf  construction-N-bigb.eps
!mv construction-N-bigb.pdf ../../LaTeX/Ola/LaTeX


set output "query-N-smallb.eps"
set nologscale y
plot "query_b=128.dat" with points pointtype 13 pointsize 2, (.6*x*x+7*x)/7000
set output
!epstopdf  query-N-smallb.eps
!mv query-N-smallb.pdf ../../LaTeX/Ola/LaTeX

set output "query-N-bigb.eps"
plot "query_b=32k.dat" with points pointtype 13 pointsize 2,  (0.007*x**2+0.08*x)
set output
!epstopdf  query-N-bigb.eps
!mv query-N-bigb.pdf ../../LaTeX/Ola/LaTeX


set output "external-constr-vs-b.eps"
set logscale x
set xlabel "b"
set nologscale y
set yrange [80:140]
plot "external_constr_N=2.dat" with points pointtype 13 pointsize 2 
set output
!epstopdf  external-constr-vs-b.eps
!mv external-constr-vs-b.pdf ../../LaTeX/Ola/LaTeX

set output "external-query-vs-b.eps"
set logscale xy
set yrange [*:*]
plot "external_query_N=2.dat" with points pointtype 13 pointsize 2 
set output
!epstopdf  external-query-vs-b.eps
!mv external-query-vs-b.pdf ../../LaTeX/Ola/LaTeX

set nologscale xy
set output "update-N-smallb.eps"
set xlabel "N"
plot "update_b=128.dat" with points pointtype 13 pointsize 2
set output
!epstopdf  update-N-smallb.eps
!mv update-N-smallb.pdf ../../LaTeX/Ola/LaTeX

set logscale xy
set output "update-b.eps"
set xlabel "b"
plot "update_N=2.dat" with points pointtype 13 pointsize 2
set output
!epstopdf  update-b.eps
!mv update-b.pdf ../../LaTeX/Ola/LaTeX

set output "time-varying-range32k.eps"
set nologscale xy
set yrange[0.23:0.25]
set xtics ("0" 0,"2e8" 2e8,"4e8" 4e8, "6e8" 6e8,"8e8" 8e8,"1e9" 1e9)
plot "rangesums_N=2_b=32k.dat" with points,\
 "firstmoments_N=2_b=32k.dat" with points
set output
!epstopdf  time-varying-range32k.eps
!mv time-varying-range32k.pdf ../../LaTeX/Ola/LaTeX

set output "time-varying-range128.eps"
set yrange [0.003:0.004]
plot "rangesums_N=2_b=128.dat" with points,\
 "firstmoments_N=2_b=128.dat" with points
set output
!epstopdf  time-varying-range128.eps
!mv time-varying-range128.pdf ../../LaTeX/Ola/LaTeX

