#!/bin/bash

path=../XAS
imp=../Input_Block.in
dollar='$'

elem=O
w=2
chg_approx=FCH
fana=../ana.out
sigma=0.4

cori_option=''
nnodes=12

if [[ $HOST == *"cori"* ]]; then
    cori_option='#SBATCH -C haswell'
fi

# get energy shifts from ana.out and import them into an array called es
declare $(grep "For element $elem" $fana | awk '{print $3, $(NF-1)}' | sed "s/^$elem/es[/" | sed 's/ /]=/')

molname=vo2
xyzname=vo2

count=0
for d in $path/$xyzname/$elem*
do
    atom=$(basename $d)

    mkdir -p $atom
    cd $atom

    index=$(echo $atom | sed "s/$elem//")
    indexw=$(seq -f "%0${w}g" $index $index)

    cat > mbxaspy.in << EOF
# initial(ground)-state
path_i          =       '../$path/$xyzname/GS'    # path
mol_name_i      =       '$molname'                    # file prefix
nbnd_i          =       800                            # number of orbitals

path_f          =       '../$path/$xyzname/$atom'    # path
mol_name_f      =       '$molname.${elem}${indexw}-$chg_approx'                    # file prefix
nbnd_f          =       800                                # number of orbitals
xas_arg         =       5                               # number of k points along one direction
gamma_only      = True
nproc_per_pool  =       1                       # processors per pool
final_1p                = True                                                  # Need one-body final-state spectrum
#spec0_only              = True
xi_analysis             = True                                                 # print out an analysis of the xi matrix
#do_paw_correction       = True

maxfn                   = 1
I_thr                   = 1e-3

ELOW                    = -5
EHIGH                   = 35
SIGMA                   = $sigma
NENER                   = 1000
ESHIFT_FINAL    = ${es[$index]}                                               # *** Please note that this is the final ESHIFT
EOF

    cat > mbxaspy.sh << EOF
#!/bin/bash
#SBATCH -p regular
#SBATCH -t 03:00:00
# SBATCH -p debug
# SBATCH -t 00:30:00
#SBATCH -e mbpy.err
#SBATCH -o mbpy.out
#SBATCH -J mbpy
#SBATCH -N $nnodes
$cori_option

MBXASPY=~/mbxaspy/main.py
srun -n 125 -c 2 --cpu_bind=cores python ${dollar}MBXASPY mbxaspy.in > mbxaspy.out

EOF

    cd ../   

    count=$((count+1))
done

cat > mbxaspy.sh << EOF
#!/bin/bash
# SBATCH -p regular
# SBATCH -t 00:10:00
#SBATCH -p debug
#SBATCH -t 00:30:00
#SBATCH -e mbpy.err
#SBATCH -o mbpy.out
#SBATCH -J mbpy
#SBATCH -N $(($nnodes*$count))
$cori_option

MBXASPY=~/mbxaspy/main.py

for atom in $elem*
do
    cd ${dollar}atom
    srun -N $nnodes -n 125 -c 2 --cpu_bind=cores python ${dollar}MBXASPY mbxaspy.in > mbxaspy.out
    cd ../
done

EOF

cat > submit_all.sh << EOF
#!/bin/bash

for atom in $elem*
do
    cd ${dollar}atom
    sbatch mbxaspy.sh
    cd ../
done
EOF

chmod 750 submit_all.sh

