BoxSize=1.5
CutOff=0.6
GMX_BIN=gmx_mpi


cat << EOF > NaCl.top
#include "oplsaa.ff/forcefield.itp"
#include "oplsaa.ff/tip3p.itp"
#include "oplsaa.ff/ions.itp"

[ System ]
NaCl 

[ Molecules ]
EOF

${GMX_BIN} solvate -cs spc216.gro   -box ${BoxSize} -o Waterbox.gro -p NaCl.top

cat << EOF > Ions.mdp
; ions.mdp - used as input into grompp to generate ions.tpr
rlist = ${CutOff}
rcoulomb = ${CutOff} 
rvdw = ${CutOff}
EOF

${GMX_BIN}  grompp -f Ions.mdp -c Waterbox.gro -p NaCl.top -o Ions.tpr
echo "2" | ${GMX_BIN}  genion -s Ions.tpr -o NaCl_IntialGeo.gro -nname CL -nn 1 -pname NA -np 1 -p NaCl.top

NumAtom_NA=`cat NaCl_IntialGeo.gro | grep NA  | awk '{print $3}'`
NumAtom_CL=`cat NaCl_IntialGeo.gro | grep CL  | awk '{print $3}'`

cat << EOF > EnergyMin.mdp
integrator = steep 
nsteps = 10000
cutoff-scheme = Verlet
coulombtype = PME
rlist = ${CutOff}
rcoulomb = ${CutOff} 
rvdw = ${CutOff}
constraints = h-bonds
EOF

${GMX_BIN}  grompp -f EnergyMin.mdp -c NaCl_IntialGeo.gro -p NaCl.top -o NaCl_EnergyMin.tpr
${GMX_BIN}  mdrun -deffnm NaCl_EnergyMin -v 
rm -f step*.pdb 

cat << EOF > MD-NPT.mdp
integrator = md 
dt = 0.002
nsteps = 1000000
cutoff-scheme = Verlet
coulombtype = PME
rlist = ${CutOff}
rcoulomb = ${CutOff}
rvdw = ${CutOff}
constraints = h-bonds
tcoupl =  V-rescale
ref_t = 300
tau-t = 1.0
tc-grps = System
gen-vel = yes
gen-temp = 300
gen-seed = -1 
DispCorr = AllEnerPres
pcoupl = Parrinello-Rahman 
pcoupltype = isotropic
ref-p = 1.01325
compressibility = 4.5e-5
nstxout-compressed = 125
nstxout = 50000
nstvout = 50000
EOF

cat << EOF > plumed.dat
ene: ENERGY
vol: VOLUME
box: CELL
dist: DISTANCE ATOMS=${NumAtom_NA},${NumAtom_CL}
UPPER_WALLS ...
   ARG=dist
   AT=0.6
   KAPPA=4000.0
   LABEL=uwall
... UPPER_WALLS
PRINT ARG=ene,dist,vol,box.ax,uwall.* FILE=colvar.data STRIDE=125
EOF

${GMX_BIN}  grompp -f MD-NPT.mdp -c NaCl_EnergyMin.gro -p NaCl.top -o NaCl_NPT.tpr -maxwarn 1 
${GMX_BIN}  mdrun -deffnm NaCl_NPT -v -plumed plumed.dat 
echo "0" | ${GMX_BIN} trjconv -f NaCl_NPT.xtc -s NaCl_NPT.tpr -pbc whole -o NaCl_NPT.pbc-whole.xtc

for i in `seq 1 10`
do 
  time=`echo "${i}*200" | bc -l`
  echo "0" | ${GMX_BIN} trjconv -f NaCl_NPT.trr -s NaCl_NPT.tpr -o NaCl_StartingStructure-${i}.gro -dump ${time} 
done 











