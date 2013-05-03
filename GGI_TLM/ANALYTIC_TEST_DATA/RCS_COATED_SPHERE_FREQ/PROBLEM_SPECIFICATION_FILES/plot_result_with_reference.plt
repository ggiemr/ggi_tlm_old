#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_4_RCS_coated_sphere_rcs_freq_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "coated.rcs" u 1:4 title "RCS: GGI_RCS_sphere" w l,\
     "PROBLEM_SPECIFICATION_FILES/coated.rcs_ref" u 1:4 title "RCS: GGI_RCS_sphere reference" w p
     
#PAUSE pause -1
