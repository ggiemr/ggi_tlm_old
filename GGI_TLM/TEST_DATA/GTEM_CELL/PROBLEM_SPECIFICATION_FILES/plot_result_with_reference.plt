#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_58_GTEM_CELL_current_ref.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/GTEM_CELL.cable_current.tout_ref" u 1:3 title "GTEM_CELL cable_current: GGI_TLM Reference " w p,\
     "GTEM_CELL.cable_current.tout" u 1:3 title "GTEM_CELL cable_current: GGI_TLM" w l 
#PAUSE pause -1
