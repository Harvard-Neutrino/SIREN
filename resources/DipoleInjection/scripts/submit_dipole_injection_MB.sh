TOT_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/tot_xsec_Enu/"
DIF_XSEC_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/diff_xsec_z_Enu/"
LI_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/"
FLUX_PATH="/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/BNB_Flux_Tables/"
FLUX_FILE="${FLUX_PATH}BNB_numu_flux.txt"
EARTH_MODEL="PREM_miniboone"
MAT_MODEL="MiniBooNE"
Z="1"
N="1e5"
SEED="7" 




cat mHNL_input_float.txt | while read M; do
cat dDipole_input.txt | while read D; do
		OUTPUT="Outputs_MB/MB_mHNL_${M}_d_${D}"
		FILE="${OUTPUT}.csv"
		echo "M = ${M} GeV; D = ${D} GeV^-1"
		if test -f "$FILE"; then
				echo "$FILE exists"
				continue
		fi
		./inject_dipole_MB --output $OUTPUT --li $LI_PATH --earth-model $EARTH_MODEL --materials-model $MAT_MODEL --tot-xsec-path $TOT_XSEC_PATH --diff-xsec-path $DIF_XSEC_PATH -z $Z -n $N -m $M -d $D --seed $SEED --numu --flux-file $FLUX_FILE 
done
done
