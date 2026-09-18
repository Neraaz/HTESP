#!/bin/bash
# Writen by Niraj K. Nepal, Ph.D.
# Script to create job submission files for different calculations

which_calc=$1
parallel_command=$2
nproc=$3
npscf=$(echo $nproc | bc)
scf_command="$parallel_command -np $npscf pw.x < scf.in > scf.out"
scf_dos="$parallel_command -np $npscf pw.x < scf-1.in > scf-1.out"
vasp_command="$parallel_command -np $npscf vasp_std"
band_command="$parallel_command -np $npscf pw.x < scf-band.in > scf-band.out"
band_process="$parallel_command -np $npscf bands.x < band.in > band.out"
dos_command="$parallel_command -np $npscf pw.x < scf-dos.in > scf-dos.out"
dos_process="$parallel_command -np $npscf dos.x < dos.in > dos.out"
pdos_process="$parallel_command -np $npscf projwfc.x < pdos.in > pdos.out"
wannier_prepare="wannier90.x -pp ex"
pw2_wann="$parallel_command -np $npscf pw2wannier90.x -in pw2wan.in > pw2wan.out"
ph_command="$parallel_command -np $nproc ph.x < elph.in > elph.out"
q2r_command="$parallel_command -np $npscf q2r.x < q2r.in > q2r.out"
dynmat_command="$parallel_command -np $npscf dynmat.x < dynmat.in > dynmat.out"
matdyn_command="$parallel_command -np $npscf matdyn.x < matdyn.in > matdyn.out"
matdyndos_command="$parallel_command -np $npscf matdyn.x < matdyn-dos.in > matdyn-dos.out"
lambda_command="$parallel_command -np $npscf lambda.x < lambda.in > lambda.out"
nscf_command="$parallel_command -np $npscf pw.x < nscf.in > nscf.out"
nscf_proj_command="$parallel_command -np $npscf pw.x < nscf-proj.in > nscf-proj.out"
projwfc_command="$parallel_command -np $npscf projwfc.x -in projwfc.in > projwfc.out"
nscf_epw="$parallel_command -np $npscf pw.x < nscf_epw.in > nscf_epw.out"
epw="$parallel_command -np $nproc epw.x -npools $nproc -i epw.in > epw.out"
wan_band="wannier90.x ex"
kadd="cat scf-band-header.in *_band.kpt > scf-band.in"
ifermi_command="ifermi plot --property velocity --interpolation-factor 10 --property-colormap bwr"

if [ $which_calc == 'help' ]; then
  echo "_______________________________________________________________________________________________"
  echo " "
  echo "Usage: generate_submission_file.sh <which_calc> <parallel_command> <number_of_processors>"
  echo "Available options: which_calc = qe-elph, epw-elph, wannier_band, and vasp. parallel_command could be mpirun, ibrun, ..."
  echo " "
  echo " "
  echo "------------------------------------------------------------------------------"
  echo "Please provide your batch submission format with name 'batch.header'"
  echo "**********************************************************************"
  echo "* #!/bin/bash                                                        *"
  echo "* #SBATCH --partition=dense --ntasks=2 --cpus-per-task=24 --time=1-0 *"
  echo "* load your necessary modules                                        *"
  echo "* submission here                                                    *"                                  
  echo "* date                                                               *"
  echo "**********************************************************************"
  echo "'submission here' line is important and should be present. This line will be replaced by actual commands"
  echo "Usage: generate_submission_file.sh qe-elph mpirun 48"
  echo "This will create submission files for electron phonon calculations in QE"
  echo "It uses nproc=ntasks*cpu-per-task = 48 cores for electron phonon calculation"
  echo "while, using only half of the cores for other small scf, nscf calculations"
elif [ $which_calc == 'qe-elph' ]; then
  sed "s/submission here/$scf_command/g" batch.header > run.sh 
  sed "s/submission here/$band_command/g" batch.header > run-band.sh 
  sed "s/submission here/$band_process/g" batch.header > run-bandp.sh 
  sed "s/submission here/$scf_dos/g" batch.header > run-dos-temp.sh 
  sed "/$scf_dos/a $dos_command" run-dos-temp.sh > run-dos.sh 
  sed "s/submission here/$dos_process/g" batch.header > run-dosp.sh 
  sed "s/submission here/$ph_command/g" batch.header > run-elph.sh 
  sed "s/submission here/$q2r_command/g" batch.header > q2r.sh 
  sed "s/submission here/$dynmat_command/g" batch.header > dynmat.sh 
  sed "s/submission here/$matdyn_command/g" batch.header > matdyn.sh 
  sed "s/submission here/$matdyndos_command/g" batch.header > matdyn-dos.sh 
  sed "s/submission here/$lambda_command/g" batch.header > lambda.sh 
  sed "s/submission here/$pdos_process/g" batch.header > pdos.sh 
  rm run-dos-temp.sh
elif [ $which_calc == 'epw-elph' ]; then
  sed "s/submission here/$scf_command/g" batch.header > temp.sh 
  sed "/$scf_command/a $ph_command" temp.sh > run-ph.sh
  sed "/$scf_command/a $nscf_proj_command" temp.sh > temp2.sh
  sed "/$nscf_proj_command/a $projwfc_command" temp2.sh > run-proj.sh
  sed "s/submission here/$nscf_epw/g" batch.header > temp3.sh
  sed "/$nscf_epw/a $epw" temp3.sh > run-epw.sh
  rm temp.sh temp2.sh temp3.sh
  # now need to add final Tc calculation


elif [ $which_calc == 'wannier_band' ]; then
  sed "s/submission here/$scf_command/g" batch.header > wannier_band.sh
  sed -i "/$scf_command/a $nscf_command" wannier_band.sh
  sed -i "/$nscf_command/a $wannier_prepare" wannier_band.sh
  sed -i "/$wannier_prepare/a $pw2_wann" wannier_band.sh
  sed -i "/$pw2_wann/a $wan_band" wannier_band.sh
  sed "s/submission here/$scf_command/g" batch.header > wannier_band2.sh
  sed -i "/$scf_command/a $band_command" wannier_band2.sh
  sed -i "/$band_command/a $band_process" wannier_band2.sh

elif [ $which_calc == 'vasp' ]; then
  sed "s/submission here/$vasp_command/g" batch.header > run_vasp.sh
  sed "s/submission here/$ifermi_command/g" batch.header > run_ifermi.sh
fi  
