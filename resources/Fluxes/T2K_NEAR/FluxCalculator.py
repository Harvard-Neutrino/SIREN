import os


def format_scientific(value):
    # Convert to float and then to scientific notation with 17 digits
    return f"{2e-16*float(value):.17e}"

def process_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        # Write header
        outfile.write("Average Energy [GeV]    numu/m^2/GeV/POT\n")

        # Skip the header lines
        for _ in range(3):
            next(infile)

        for line in infile:
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                energy_range = parts[1:4]
                numu = parts[4]

                # Calculate average energy
                start_energy = float(energy_range[0])
                end_energy = float(energy_range[2])
                avg_energy = (start_energy + end_energy) / 2

                # Format both average energy and numu to scientific notation with 17 digits
                formatted_avg_energy = format_scientific(avg_energy)
                formatted_numu = format_scientific(numu)

                # Write to .dat file
                outfile.write(f"{formatted_avg_energy}    {formatted_numu}\n")

# Usage
input_filename = 't2kflux_2020_plus250kA_runcond_nd280.txt'  # Your input .txt filename
output_filename = 'T2KFLUX.dat'  # Your desired output .dat filename

process_file(input_filename, output_filename)
print(f"Data has been processed and saved to {output_filename}")