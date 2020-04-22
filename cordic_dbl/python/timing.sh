#! /bin/bash

# ./timing.sh | tee runtime.txt | awk '/i_Ecs/{printf "%-10g", $2} {printf " %-6g",$9=="usec"?$8/1000.:$8}  /ke._E_newton_my/{print ""}' | awk '{getline x; print $0,x}' | tee runtime.dat


for e in 0. `seq 0 0.05 0.95` 0.999999 ; do
   echo $e >& 2
   for setup in "M=np.arange(0, np.pi, np.pi/10000);" "E = np.arange(0,np.pi,np.pi/10000); M = E - $e*np.sin(E);"; do
      printf "%-17s %s " "i_Ecs" $e      # CORDIC double rotation integer
      python -m timeit -s "import numpy as np, ke_dbl_i; $setup"  "ke_dbl_i.Ecs(M, "$e", typ='i_Ecs')"

      printf "%-17s %s " "f_Ecs" $e      # CORDIC double rotation float
      python -m timeit -s "import numpy as np, ke_dbl_f; $setup"  "ke_dbl_f.Ecs(M, "$e", typ='f_Ecs')"

      printf "%-17s %s " "ke._E" $e      # CORDIC like
      python -m timeit -s "import numpy as np, ke    ; $setup"  "ke._E(M, "$e", n=29)"

      printf "%-17s %s " "ke._E_newton" $e   # Newton with built-in sin
      python -m timeit -s "import numpy as np, ke    ; $setup"  "ke._E_newton(M, "$e")"

      printf "%-17s %s " "ke._E_newton_my" $e  # Newton with own sin
      python -m timeit -s "import numpy as np, ke    ; $setup"  "ke._E_newton(M, "$e", typ='my')"
   done
done


