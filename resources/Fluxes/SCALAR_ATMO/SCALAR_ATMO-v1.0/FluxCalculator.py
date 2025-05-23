import os
def MakeFluxFile(tag, abs_flux_dir):
    '''
    Accepts the following tags:
        phi
    '''
    mode,particle = tag.split("_")
    if mode not in ["phi"]:
        print("%s beam mode specified in tag %s is not valid"%(mode,tag))
        exit(0)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "%s.dat"%mode)
    output_flux_file = os.path.join(abs_flux_dir,
                                    "%s.txt"%mode)
    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        headers = all_lines[0].strip().split()
        data = [line.strip().split() for line in all_lines[1:]]
        pid = headers.index(particle)
        with open(output_flux_file,"w") as fout:
            for row in data:
                energy,flux = float(row[0]),float(row[1])
                flux = flux*1e4*4*3.14159/(energy**3) # put flux in units of nu/m^2/GeV/s
                print(energy,flux,file=fout)
    return output_flux_file