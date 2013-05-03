#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_58_GTEM_CELL_current.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"
plot "GTEM_CELL.cable_current.tout" u 1:3 title "GTEM_CELL cable_current: GGI_TLM" w l 
#PAUSE pause -1
