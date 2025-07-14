import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
"""
This script processes X-ray diffraction (XRD) data for non-uniform layers, 
applying corrections to the layer mixture, calculating the index of refraction, 
attenuation coefficients, and corrected depths. It reads data from a excel file 
previously created, performs calculations, and stores the results in a new excel file.

Mark Smith
"""

def get_data(file_path: str):
    """
    Reads data from a excel file and returns it as a dictionary.
    Args:
        file_path (str): The path to the excel file.
    Returns:
        dict: The data read from the excel file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    data_frame = pd.read_excel(file_path)

    data = {'file_name': data_frame['Unnamed: 0'],
            'omega': data_frame['omega'].to_numpy(),
            'amorphous_percentage': data_frame['amorphous_percentage'].to_numpy(),
            'crystalline_percentage': data_frame['crystalline_percentage'].to_numpy(),
            'effective_depth': data_frame['effective_depth'].to_numpy() * 10**-6,
            'amorphous_percentage_intensity': data_frame['amorphous_percentage_intensity'].to_numpy(),
            'crystalline_percentage_intensity': data_frame['crystalline_percentage_intensity'].to_numpy(),
            'effective_depth_intensity': data_frame['effective_depth_intensity'].to_numpy() * 10**-6}

    return data

def Layer_Mixture_Correction(data: dict):
    """
    Applies a correction to the layer mixture data.
    Args:
        data (dict): The data dictionary containing layer mixture information.
    Returns:
        dict: The corrected data dictionary.
    """
    # Unpack the data
    amorphous_percentage = data['amorphous_percentage']
    crystalline_percentage = data['crystalline_percentage']
    effective_depth = data['effective_depth']

    amorphous_percentage_intensity = data['amorphous_percentage_intensity']
    crystalline_percentage_intensity = data['crystalline_percentage_intensity']
    effective_depth_intensity = data['effective_depth_intensity']

    # Weighting factors
    weight = [(effective_depth[i] / effective_depth[i+1]) for i in range(len(effective_depth)-1)]
    
    weight_intensity = [(effective_depth_intensity[i] / effective_depth_intensity[i+1]) for i in range(len(effective_depth_intensity)-1)]
    
    # Calculate the corrected mixture values
    corrected_amorphous_percentage = [(amorphous_percentage[i+1] - (amorphous_percentage[i] * weight[i])) / (1 - weight[i]) for i in range(len(amorphous_percentage)-1)]
    corrected_crystalline_percentage = [(crystalline_percentage[i+1] - (crystalline_percentage[i] * weight[i])) / (1 - weight[i]) for i in range(len(crystalline_percentage)-1)]

    corrected_amorphous_percentage_intensity = [(amorphous_percentage_intensity[i+1] - (amorphous_percentage_intensity[i] * weight_intensity[i])) / (1 - weight_intensity[i]) for i in range(len(amorphous_percentage_intensity)-1)]
    corrected_crystalline_percentage_intensity = [(crystalline_percentage_intensity[i+1] - (crystalline_percentage_intensity[i] * weight_intensity[i])) / (1 - weight_intensity[i]) for i in range(len(crystalline_percentage_intensity)-1)]

    # Add the first element back to the corrected lists
    corrected_amorphous_percentage.insert(0, amorphous_percentage[0])
    corrected_crystalline_percentage.insert(0, crystalline_percentage[0])
    
    corrected_amorphous_percentage_intensity.insert(0, amorphous_percentage_intensity[0])
    corrected_crystalline_percentage_intensity.insert(0, crystalline_percentage_intensity[0])

    # Add in a first value for the weights, so that the lengths match
    weight.insert(0, 1.0)
    weight_intensity.insert(0, 1.0)

    # Update the data dictionary with corrected values
    data['corrected_amorphous_percentage'] = corrected_amorphous_percentage
    data['corrected_crystalline_percentage'] = corrected_crystalline_percentage
    data['weight'] = weight

    data['corrected_amorphous_percentage_intensity'] = corrected_amorphous_percentage_intensity
    data['corrected_crystalline_percentage_intensity'] = corrected_crystalline_percentage_intensity
    data['weight_intensity'] = weight_intensity
    
    
    return data

def Index_Of_Refraction(data: dict):
    """
    Calculates the index of refraction based on the material type and the corrected layer mixture data.
    Args:
        data (dict): The data dictionary containing layer mixture information.
    Returns:
        dict: The data dictionary with the calculated index of refraction.
    """

    # Unpack the data
    mixture_amorphous = data['corrected_amorphous_percentage']
    mixture_crystalline = data['corrected_crystalline_percentage']

    mixture_amorphous_intensity = data['corrected_amorphous_percentage_intensity']
    mixture_crystalline_intensity = data['corrected_crystalline_percentage_intensity']

    file_name = data['file_name'][0] # Only one file name is needed for determining the material type

    # website for index of refraction https://henke.lbl.gov/optical_constants/getdb2.html
    # Get indices of refraction based on the file name
    material = input("Enter the material type (e.g., 'P' for Phosphorus, 'B' for Boron. Case sensitive): ").strip()  # Ask the user for the material type
    if material == 'P':
        n = 1-6.96137067E-06-1.98567363E-07J # Phosphorus index of refraction
    elif material == 'B':
        n = 1-6.94444225E-06-5.98799899E-09J # Boron index of refraction
    else:
        print("Please update code to include the correct indices of refraction for your material.")
        raise ValueError("Material index of refraction unknown.")
    
    Sin = 1-7.57442876E-06-1.72761077E-07J

    # Calculate the corrected indices of refraction
    index_of_refraction = []
    index_of_refraction_intensity = []
    for index in range(len(mixture_amorphous)):
        index_of_refraction.append(mixture_amorphous[index] * n + mixture_crystalline[index] * Sin)
        index_of_refraction_intensity.append(mixture_amorphous_intensity[index] * n + mixture_crystalline_intensity[index] * Sin)
    
    # Add the indices of refraction to the data dictionary
    data['index_of_refraction'] = index_of_refraction
    data['index_of_refraction_intensity'] = index_of_refraction_intensity

    return data

def Attenuation_Coefficient(data: dict):
    """
    Calculates the attenuation coefficient based on the index of refraction.
    Args:
        data (dict): The data dictionary containing index of refraction information.
    Returns:
        dict: The data dictionary with the calculated attenuation coefficients.
    """
    
    def attenuation_coefficient(ni: complex): # Attenuation coefficient calculation (only used in this function)
        """
        Calculate the attenuation coefficient for a given material.

        Parameters:
            ni (complex): Complex refractive index of the material.

        Returns:
            float: Attenuation coefficient in m^-1.
        """
        λ = 1.5406e-10  # Cu Kα radiation wavelength in meters
        return 4 * np.pi * -ni.imag / λ  # Attenuation coefficient in m^-1

    # Unpack the data
    index_of_refraction = data['index_of_refraction']
    index_of_refraction_intensity = data['index_of_refraction_intensity']

    attenuation_coefficient_list = [] # List to store skin depth values
    attenuation_coefficient_list_intensity = [] # List to store skin depth values for intensity

    for current_index_of_refraction_index in range(len(index_of_refraction)):
        
        # Calculate the attenuation coefficient
        effective_alpha = (attenuation_coefficient(index_of_refraction[current_index_of_refraction_index]))
        effective_alpha_intensity = (attenuation_coefficient(index_of_refraction_intensity[current_index_of_refraction_index]))

        # Store them
        attenuation_coefficient_list.append(effective_alpha)        
        attenuation_coefficient_list_intensity.append(effective_alpha_intensity)

    # Store the skin depths in the data dictionary
    data['attenuation_coefficient'] = attenuation_coefficient_list
    data['attenuation_coefficient_intensity'] = attenuation_coefficient_list_intensity

    return data

def Depth_Calculation(data: dict):
    """
    Calculates the corrected depth for each mixture based on the attenuation coefficient and angle.
    Args:
        data (dict): The data dictionary containing attenuation coefficients and angles.
    Returns:
        dict: The data dictionary with the corrected depth values added.
    """

    def Irradiance(alpha: float, distance: float):
        """
        Calculate the irradiance at a given distance based on the attenuation coefficient.
        Parameters:
            alpha (float): Attenuation coefficient in m^-1.
            distance (float): Distance traveled in the medium in meters.
        Returns:
            float: Irradiance at the given distance.
        """

        irradiance_at_depth = np.exp(-alpha * distance)

        return irradiance_at_depth

    def Loop_For_Calculation(omega_list: list, depths: list, attenuation_coefficient: list):
        """
        Loops through the omega list and calculates the corrected depth for each angle.
        Parameters:
            omega_list (list): List of angles in degrees.
            depths (list): List of effective depths in micrometers.
            attenuation_coefficient (list): List of attenuation coefficients in m^-1.
        Returns:
            list: List of corrected depths in micrometers.
        """
        
        corrected_depth = []

        for omega in omega_list:
            print(f"Processing {omega} file") # Assure the user it is working
            # Setup for loop
            path_length_check = np.array(depths) / np.sin(np.radians(omega)) # This is to be able to check how far into the sample I am so I can know which attenuation coefficient to use
            index = 0 # Index to switch attenuation coefficients, and path length check
            current_path_length_check_value = path_length_check[index]
            dy = 1e-10 # Infinitesimal distance moved into the sample by the x-ray (in SI units)
            intensity = 1 # Starting intensity is 1 (100%)
            alpha = attenuation_coefficient[index] # Starting attenuation coefficient
            depth = dy # Starting depth of x-ray
            intensity_attenuation = Irradiance(alpha, dy) # Calculate starting intensity attenuation

            while intensity > 0.1: # Go until 90% attenuation

                intensity *= intensity_attenuation
                depth += dy # Update to a new depth

                if depth >= current_path_length_check_value: # Check to see if we have gone down to a new layer
                    index += 1 # Update index

                    if index > len(path_length_check): # Stop if index exceeds length
                        break

                    current_path_length_check_value = path_length_check[index] # Update new check requirment
                    alpha = attenuation_coefficient[index] # Update new alpha value
                    intensity_attenuation = Irradiance(alpha, dy) # Calculate new intensity attenuation

            corrected_depth.append(depth * np.sin(np.radians(omega)) * 1e6) # Convert to micrometers

        return corrected_depth

    # Unpack data
    attenuation_coefficient = data['attenuation_coefficient']
    attenuation_coefficient_intensity = data['attenuation_coefficient_intensity']

    depths = data['effective_depth']
    depths_intensity = data['effective_depth_intensity']

    omega_list = data['omega']
    
    corrected_depth = Loop_For_Calculation(omega_list, depths, attenuation_coefficient)
    corrected_depth_intensity = Loop_For_Calculation(omega_list, depths_intensity, attenuation_coefficient_intensity)
    
    data['corrected_depth'] = corrected_depth
    data['corrected_depth_intensity'] = corrected_depth_intensity

    return data

def Store_Data(data: dict, path_name: str):
    """
    Stores the processed data into a excel file.
    
    Parameters:
        data (dict): A dictionary containing processed XRD data.
        path_name (str): The path where the excel file will be saved.
    """
    try:
        df = pd.DataFrame.from_dict(data)  # Convert the dictionary to a DataFrame
        df = df.sort_values(by='omega', ascending=True)  # Sort by Omega from smallest to largest
        df.to_excel(os.path.join(path_name, 'Corrected_Depth_XRD_Data.xlsx'))  # Save the DataFrame to a excel file
    
    except PermissionError as e:
        print(f"Permission denied: {e}. Please check if the file is open, and close it.")
        input("Press Enter to try saving again after closing the file.")
        df.to_excel(os.path.join(path_name, 'Corrected_Depth_XRD_Data.xlsx'))  # Save the DataFrame to a excel file

def Corrected_Depths_For_Nonuniform_Layers(path_name: str, path_name_for_storage: str):
    """
    Reads data from a excel file and applies corrections for non-uniform layers.
    Parameters:
        path_name (str): The path to the excel file.
    """
    
    data = get_data(path_name)
    
    corrected_mixture = Layer_Mixture_Correction(data)

    index_of_refraction = Index_Of_Refraction(corrected_mixture)
    
    attenuation_coefficient_data = Attenuation_Coefficient(index_of_refraction)

    finished_data = Depth_Calculation(attenuation_coefficient_data)

    Store_Data(finished_data, path_name_for_storage)

def bar_graph_depths(path_name: str, area: str = ''):
    """
    Create a bar graph of the depths calculated for each mixture.
    """

    df = pd.read_excel(path_name + '//Corrected_Depth_XRD_Data.xlsx')
    x = np.arange(len(df['omega']))
    plt.bar(x, df['corrected_crystalline_percentage'+ area] + df['corrected_amorphous_percentage'+ area], label="Si")
    plt.bar(x, df['corrected_amorphous_percentage'+ area], label="P")
    plt.xlabel('Omega (°)')
    plt.ylabel('Mix')
    plt.title('Depths for Mixtures'+ area)
    plt.xticks(ticks=x, labels=[f"{df['omega'][i]}ω-{df['corrected_depth'+ area][i]:.2f}μm"for i in range(len(df['omega']))], rotation=90)
    plt.legend()
    plt.ylim(-.05, 1.05)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    material = input("p, p01, or dt: ")
    storage_path = f'E://Full//{material} csv//' # This should be the folder that has the excel to be corrected
    # storage_path = input("Enter the path to the XRD data folder: ").strip()  # Get the path to the XRD data files from the user
    file_path = storage_path + 'Processed_XRD_Data.xlsx' # This is the name of the excel that will be corrected (this is the default name)

    try:
        Corrected_Depths_For_Nonuniform_Layers(file_path, storage_path)
        bar_graph_depths(storage_path, area = '_intensity') # Create a bar graph of the depths calculated for each mixture
        # bar_graph_depths(storage_path, area = '_sum') # Create a bar graph of the depths calculated for each mixture
        bar_graph_depths(storage_path) # Create a bar graph of the depths calculated for each mixture
    except FileNotFoundError as e:
        print(e)