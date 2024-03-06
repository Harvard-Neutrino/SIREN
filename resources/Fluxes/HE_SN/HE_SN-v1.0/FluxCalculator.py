import os

def MakeFluxFile(tag, abs_flux_dir):
    '''
    only supported tag is "numu"
    '''
    if tag!="numu":
        print("Tag %s not supported for HE SN"%tag)
        exit(0)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "dN_dE_SNe_2n_D1_0_s20_t100d_NuMu_d10kpc.txt")
    output_flux_file = os.path.join(abs_flux_dir,
                                    "HE_SN_numu.txt")
    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        data = [line.strip().split() for line in all_lines]
        with open(output_flux_file,"w") as fout:
            for row in data:
                E,flux = float(row[0]),float(row[1])
                flux*=1e4 # put flux in units of nu/m^2/GeV/100d
                print(E,flux,file=fout)
    return output_flux_file 