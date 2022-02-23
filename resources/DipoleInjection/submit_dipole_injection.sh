TOT_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/"
DIF_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/diff_xsec_z_Enu/"
LI_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/"
EARTH_MODEL="PREM_minerva"
MAT_MODEL="Minerva"
Z="1"
N="1e3"
HNL_MASSES={"0.1","0.4"} #GeV
D_DIPOLES={"1e-7","3e-7","5e-7"} #GeV^-1
SEED="7" 


# LE LOOP
for beam in {"FHC","RHC"}; do
for primary in {"nue","nuebar","numu","numubar"}; do
for M in {"0.1","0.4"}; do
for D in {"1e-7","3e-7","5e-7"}; do
		OUTPUT="LE_${beam}_${primary}_mHNL_${M}_d_${D}"
		echo $OUTPUT
		#echo ./inject_dipole --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --${beam} --LE --${primary} --seed $SEED 
		./inject_dipole --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --${beam} --LE --${primary} --seed $SEED 
done
done
done
done


# ME FHC 
beam="FHC"
primary="numu"
for M in {"0.1","0.4"}; do
for D in {"1e-7","3e-7","5e-7"}; do
		OUTPUT="ME_${beam}_${primary}_mHNL_${M}_d_${D}"
		echo $OUTPUT
		./inject_dipole --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --${beam} --ME --${primary} --seed $SEED 
done
done
