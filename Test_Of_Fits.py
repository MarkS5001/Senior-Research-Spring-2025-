import os
import pandas as pd
import matplotlib.pyplot as plt
from Peak_Area_Fitting import Gaussian_function
import numpy as np

def get_data(file_path: str):
    """
    Reads data from a CSV file and returns it as a dictionary.
    Args:
        file_path (str): The path to the CSV file.
    Returns:
        dict: The data read from the CSV file.
    """
    data = {}

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    data_frame = pd.read_excel(file_path)

    return data_frame

def Fit_Test(data: pd.DataFrame):

    for index, row in data.iterrows():
        try:
            # Unpack data for each file
            position = np.fromstring(row['Position'].strip('[]'), sep=' ') # Convert string representation of list to numpy array
            intensity = np.fromstring(row['Intensity'].strip('[]'), sep=' ') # Convert string representation of list to numpy array
            amorphous_amplitude = row['amorphous_fit_amplitude']
            amorphous_mean = row['amorphous_fit_mean']
            amorphous_std_dev = row['amorphous_fit_std_dev']
            amorphous_background = row['amorphous_fit_background']
            # amorphous_x_range = row['amorphous_fit_amorphous_x_range_used']
            crystalline_amplitude = row['crystalline_fit_amplitude']
            crystalline_mean = row['crystalline_fit_mean']
            crystalline_std_dev = row['crystalline_fit_std_dev']
            crystalline_background = row['crystalline_fit_background']
            # crystalline_x_range = row['crystalline_fit_crystalline_x_range_used']

            if len(position) != len(intensity):
                if len(intensity) > len(position):
                    intensity = intensity[:len(position)]  # Ensure intensity matches position length
                elif len(position) > len(intensity):
                    position = position[:len(intensity)]

            fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True) # Make plots to compare, and share y-axis for easier comparison
            # Plot it
            ax1.plot(position, intensity, label='Intensity', color='black')
            ax1.plot(position, Gaussian_function(position, amorphous_amplitude, amorphous_mean, amorphous_std_dev, amorphous_background), label='Amorphous Fit', color='blue')
            ax1.plot(position, Gaussian_function(position, crystalline_amplitude, crystalline_mean, crystalline_std_dev, crystalline_background), label='Crystalline Fit', color='red')

            ax2.plot(position, intensity, label='Intensity', color='black')
            ax2.plot(position, Gaussian_function(position, amorphous_amplitude, amorphous_mean, amorphous_std_dev), label='Amorphous Fit', color='blue')
            ax2.plot(position, Gaussian_function(position, crystalline_amplitude, crystalline_mean, crystalline_std_dev), label='Crystalline Fit', color='red')

            # Make it pretty
            fig.supxlabel('Position') # Set x-label for the whole figure
            fig.supylabel('Intensity') # Set y-label for the whole figure
            fig.suptitle(row['Unnamed: 0']) # Set a single title for the whole figure
            ax2.set_yscale('log') # Set log scale

            plt.tight_layout()
            plt.legend()
            plt.show()
        except Exception as e:
            print(row['Unnamed: 0'])
            print(e)
            plt.close(fig)

if __name__ == "__main__":

    # Get the file path
    material = input("p or dt: ")
    storage_path = f'E://Full//{material} csv//' # This should be the folder that has the file
    file_path = storage_path + 'Processed_XRD_Data.xlsx' # This is the name of the csv that will be corrected (this is the default name)
    
    try:
        # Get the data
        data= get_data(file_path)
        Fit_Test(data)        

    except FileNotFoundError as e: # Error handling for file not found
        print(e)