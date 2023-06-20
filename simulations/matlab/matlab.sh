#!/bin/bash
exec=${INC}
module load matlab/r2017b

matlab_script=Multi_Cracks_23_6_11_8

#python ${exec}/indent_xy.py
#-nojvm -nosplash -nodesktop -r
matlab -nodisplay -r "try, run('${exec}/${matlab_script}'), catch e, disp(getReport(e)), exit(1), end, exit(0);"
echo "matlab exit code: $?"

