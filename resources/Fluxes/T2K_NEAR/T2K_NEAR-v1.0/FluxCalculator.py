import pandas as pd
import os

def MakeFluxFile(tag, abs_flux_dir):
    # Read the file, skipping the header rows
    
    input_file=os.path.join(abs_flux_dir,'t2kflux_2020_plus250kA_runcond_nd280.txt')
    
    df = pd.read_csv(input_file, sep='\s+', skiprows=3)

    # Process the 'Energy [GeV]' column
    df['Energy_Start'] = df['Energy'].str.split('-').str[0].astype(float)
    df['Energy_End'] = df['Energy'].str.split('-').str[1].astype(float)
    df['Energy_Avg'] = (df['Energy_Start'] + df['Energy_End']) / 2

    # Process the 'numu' column
    df['numu_processed'] = df['numu'] * 2e-16

    # Create the output dataframe
    output_df = pd.DataFrame({
        'Energy_Avg': df['Energy_Avg'].round(8),
        'numu_processed': df['numu_processed'].round(8)
    })

    # Define the output file path
    output_file = os.path.join(abs_flux_dir,"t2k_flux.txt")

    # Write the output file
    output_df.to_csv(output_file, sep=' ', index=False, header=False, float_format='%.8f')

    return output_file