import os
def load_flux(tag):
    '''
    Accepts the following tags:
        {model}_{osc_status}_{particle1}_{particle2}_...{particleN}
        where model is one of {bartol,daemonflux,H3a_SIBYLL21,H3a_SIBYLL23C,honda2006}
        osc_status is one of {osc,unosc}
        particlei is one of {nue,numu,nuebar,numubar,nutau,nutaubar}
    '''
    fields = tag.split("_")
    if len(fields) < 3:
        print("Tag %s is not valid. It should be of the form {model}_{osc_status}_{particle1}_{particle2}_...{particleN}"%tag)
        exit(0)
    model, osc_status = fields[0], fields[1]
    particles = fields[2:]
    if model not in ["bartol","daemonflux","H3a_SIBYLL21","H3a_SIBYLL23C","honda2006"]:
        print("%s model specified in tag %s is not valid"%(model,tag))
        exit(0)
    if osc_status not in ["osc","unosc"]:
        print("%s osc_status specified in tag %s is not valid"%(osc_status,tag))
        exit(0)
    if not all(p in ["nue","numu","nuebar","numubar","nutau","nutaubar"] for p in particles):
        print("%s particles specified in tag %s are not valid"%(particles,tag))
        exit(0)

    abs_flux_dir = os.path.dirname(__file__)
    output_flux_file = os.path.join(abs_flux_dir,
                                    "Atmospheric_%s_flux.txt"%(tag))
    if os.path.isfile(output_flux_file):
        return output_flux_file
    else:
        with open(output_flux_file,"w") as fout:
            tot_flux = {}
            key_list = []
            for particle in particles:
                input_flux_file = os.path.join(abs_flux_dir,
                                               "%s/%s_%s_flux.txt"%(model,particle,osc_status))
                with open(input_flux_file,"r") as fin:
                    all_lines = fin.readlines()
                    data = [line.strip().split() for line in all_lines]

                    for row in data:
                        Energy,CosTheta,flux = float(row[0]),float(row[1]),float(row[2])
                        Energy_key = "%3.3e"%Energy
                        CosTheta_key = "%3.3e"%CosTheta
                        key = (Energy_key,CosTheta_key)
                        if key not in tot_flux:
                            tot_flux[key] = 0.
                            key_list.append(key)
                        tot_flux[key] += flux
            for key in key_list:
                Energy_key, CosTheta_key = key
                Energy = float(Energy_key)
                CosTheta = float(CosTheta_key)
                flux = tot_flux[key] # already in nu/m^2/GeV/POT
                print(Energy,CosTheta,flux,file=fout)
    return output_flux_file