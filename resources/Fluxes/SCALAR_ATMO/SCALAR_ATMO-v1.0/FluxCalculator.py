import os
s_per_year = 60*60*24*365.25
def MakeFluxFile(tag, abs_flux_dir):
    '''
    Accepts the following tags:
        phi
    '''
    if tag not in ["phi"]:
        print(" tag %s is not valid"%(tag))
        exit(0)
    input_flux_file = os.path.join(abs_flux_dir,
                                   "%s.dat"%tag)
    output_flux_file = os.path.join(abs_flux_dir,
                                    "%s.txt"%tag)
    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        headers = all_lines[0].strip().split()
        data = [line.strip().split() for line in all_lines[1:]]
        with open(output_flux_file,"w") as fout:
            for row in data:
                energy,flux = float(row[0]),float(row[1])
                flux = flux*s_per_year*1e4*4*3.14159/(energy**3) # put flux in units of nu/m^2/GeV/yr
                print(energy,flux,file=fout)
    return output_flux_file