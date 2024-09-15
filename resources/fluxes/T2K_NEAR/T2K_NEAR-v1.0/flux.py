import os

def load_flux(tag):

    '''
    Accepts the following tags:
        {PLUS, MINUS}_{nue,nuebar,numu,numubar}
    '''

    enhance, particle = tag.split("_")

    if enhance not in ["MINUS", "PLUS"]:
        print("%s 250kA enhancement specified in tag %s is not valid"%(enhance))
        exit(0)
    if particle not in ["numu", "numubar", "nue", "nuebar"]:
        print("%s particle specified in tag %s is not valid"%(particle))
        exit(0)

    abs_flux_dir = os.path.dirname(__file__)

    input_flux_file = os.path.join(abs_flux_dir,
                                   "T2K_%s_250kA.dat"%(enhance))

    output_flux_file = os.path.join(abs_flux_dir,
                                    "T2KOUT_%s.dat"%(tag))

    with open(input_flux_file,"r") as fin:
        all_lines = fin.readlines()
        headers = all_lines[1].strip().split()
        data = [line.strip().split() for line in all_lines[3:]]
        pid = headers.index(particle)
        with open(output_flux_file,"w") as fout:
            for row in data:
                E, flux = (float(row[1])+float(row[3]))/2, float(row[pid+2])
                flux*=2e-16 # put flux in units of nu/m^2/GeV/POT
                print(E, flux, file=fout)
    return output_flux_file
