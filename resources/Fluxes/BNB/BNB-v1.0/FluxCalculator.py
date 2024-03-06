import os
def MakeFluxFile(tag, abs_flux_dir):
    '''
    Accepts the following tags:
        {FHC,RHC}_{nue,nuebar,numu,numubar}
    '''
    mode,particle = tag.split("_")
    if mode not in ["FHC","RHC"]:
        print("%s beam mode specified in tag %s is not valid"%(mode,tag))
        exit(0)
    if particle not in ["nue","numu","nuebar","numubar"]:
        print("%s particle specified in tag %s is not valid"%(particle,tag))
        exit(0)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "BNB_%s.dat"%mode)
    output_flux_file = os.path.join(abs_flux_dir,
                                    "BNB_%s_%s_flux.txt"%(mode,particle))
    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        headers = all_lines[0].strip().split()
        data = [line.strip().split() for line in all_lines[1:]]
        pid = headers.index(particle)
        with open(output_flux_file,"w") as fout:
            for row in data:
                Elow,Ehigh,bin_flux = float(row[0]),float(row[1]),float(row[pid])
                Emid = (Elow+Ehigh)/2.
                flux = bin_flux/50*1000*1e4 # put flux in units of nu/m^2/GeV/POT
                print(Emid,flux,file=fout)
    return output_flux_file 