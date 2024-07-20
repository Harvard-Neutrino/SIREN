import os

def MakeFluxFile(tag, abs_flux_dir):
    def format_scientific(value):
        # Convert to float and then to scientific notation with 17 digits
        return f"{float(value):.17e}"

    result = "Average Energy [GeV]    numu/m^2/GeV/POT\n"

    input_flux_file=os.path.join(abs_flux_dir,'t2kflux_2020_plus250kA_runcond_nd280.txt')
    output_flux_file = os.path.join(abs_flux_dir,"t2k_flux.txt")

    
    with open(input_flux_file, 'r') as infile:
        # Skip the header lines
        for _ in range(3):
            next(infile)

        for line in infile:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                energy_range = parts[1:4]
                numu = float(parts[4])  # Convert numu to float

                # Calculate average energy
                start_energy = float(energy_range[0])
                end_energy = float(energy_range[2])
                avg_energy = (start_energy + end_energy) / 2

                # Multiply numu by 2e-16
                numu *= 2e-16

                # Format both average energy and numu to scientific notation with 17 digits
                formatted_avg_energy = format_scientific(avg_energy)
                formatted_numu = format_scientific(numu)

                # Append to result string
                result += f"{formatted_avg_energy}    {formatted_numu}\n"
                
                with open('t2k_flux.txt', 'w') as outfile:
                    outfile.write(result)

    return output_flux_file

# Usage
# input_filename = 't2kflux_2020_plus250kA_runcond_nd280.txt'  # Your input .txt filename
# output_text = MakeFluxFile(input_filename)

# Print the first few lines of the output (for verification)
# print("First few lines of the output:")
# print("\n".join(output_text.split("\n")[:5]))

# If you want to save the output to a file, you can do:
# with open('output.dat', 'w') as outfile:
#     outfile.write(output_text)