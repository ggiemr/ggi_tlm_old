#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_59_HUGO_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "hugo.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1
