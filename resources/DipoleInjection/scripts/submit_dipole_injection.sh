TOT_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/tot_xsec_Enu/"
DIF_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/diff_xsec_z_Enu/"
LI_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/"
FLUX_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/NUMI_Flux_Tables/"
EARTH_MODEL="PREM_minerva"
MAT_MODEL="Minerva"
Z="1"
N="1e5"
SEED="7" 


# LE LOOP
#for beam in {"FHC","RHC"}; do
#for primary in {"nue","nuebar","numu","numubar"}; do
#for M in {"0.1","0.4"}; do
#for D in {"1e-7","3e-7","5e-7"}; do
#		OUTPUT="LE_${beam}_${primary}_mHNL_${M}_d_${D}"
#		FLUX_FILE="${FLUX_PATH}LE_${beam}_${primary}.txt"
#		echo $OUTPUT
#		./inject_dipole --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --${beam} --LE --${primary} --seed $SEED --flux-file $FLUX_FILE 
#done
#done
#done
#done


# ME FHC 
beam="FHC"
primary="numu"
cat mHNL_input_float_abr.txt | while read M; do
cat dDipole_input.txt | while read D; do
#for D in {"1e-7","3e-7","5e-7"}; do
		OUTPUT="ME_${beam}_${primary}_mHNL_${M}_d_${D}"
		FLUX_FILE="${FLUX_PATH}ME_${beam}_${primary}.txt"
		FILE="Outputs/${OUTPUT}.csv"
		echo "M = ${M} GeV; D = ${D} GeV^-1"
		if test -f "$FILE"; then
				echo "$FILE exists"
				continue
		fi
		./inject_dipole --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --${beam} --ME --${primary} --seed $SEED --flux-file $FLUX_FILE 
done
done
