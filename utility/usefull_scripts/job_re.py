import os
import sys


celli = int(int(sys.argv[1]) + 6)
fc = int(sys.argv[2])
cellj = int(celli + 1)
os.system("""sed -n '/&CONTROL/,/CELL_PARAMETERS/p' R{}/scf.in | sed '$d' > R{}/scf-header.in""".format(fc,fc))
os.system("""cp R{}/scf.in R{}/scf-norelax.in""".format(fc,fc))
os.system("""cp R{}/scf.out R{}/scf-norelax.out""".format(fc,fc))
os.system("""grep -A {} 'CELL_PARAMETERS (angstrom)' R{}/scf.out | tail -n {} | sed '$d' > R{}/cell_norelax""".format(celli,fc,cellj,fc))
os.system("""cat R{}/scf-header.in R{}/cell_norelax > R{}/scf.in""".format(fc,fc,fc))
os.system("cp run.sh R{}/run-{}-re.sh".format(fc,fc))
os.chdir("R{}".format(fc))
os.system("sbatch run-{}-re.sh".format(fc))
os.system("sleep 2s")
os.chdir("../")
