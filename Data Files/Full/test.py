# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def XRD_Data_Dictionary(path_name: str): # Gets position and intensity from file to a dictionary
    """
    Reads XRD data from a specified path and returns a dictionary with the data.
    
    Parameters:
    path_name (str): The path to the directory containing the XRD data files.
    
    Returns:
    dict: A dictionary where keys are file names, and values are position and intensity.
    """
    xrd_data = {} # Initialize an empty dictionary to store the data
    
    os.chdir(path_name) # Change the current working directory to the specified path

    for file in os.listdir(path_name): # Iterate through each file in the directory
        if file.endswith('.csv') and 'MSS' in file: # Check if the file is a CSV file, and is a scan file
            
            scan = pd.read_csv(file, index_col=0) # Read the CSV file into a DataFrame

            intensity = scan['Intensity [Counts]'].to_numpy()  # Get the intensity data from the CSV file
            position = scan['Pos. [°2θ]'].to_numpy()  # Get the position data from the CSV file

            # Extract the omega value from the file name
            index_location_start = file.index(',')
            index_location_end = file.index('w')
            omega = float(file[index_location_start+1:index_location_end].strip())  # Extract the omega value from the file name

            xrd_data[file] = {'omega': omega, 'Position': position, 'Intensity': intensity}  # Store the data in the dictionary

    return xrd_data

def Length_Peak(data: dict): # Just gets the peak ranges
    """
    Prompts the user to input the amorphous and crystalline ranges for each file in 
    the data dictionary.
    
    Parameters:
    data (dict): A dictionary containing XRD data with file names as keys and 
    position/intensity as values.
    
    Returns:
    dict: The updated dictionary with amorphous and crystalline ranges added for each file.
    """
    # Change this later
    manual_input = 'no' # input("Do you want to manually input the ranges? (yes/no): ").strip().lower()  # Ask the user if they want to manually input the ranges
    if manual_input == 'no':
        peak_range_file = 'Peak_Ranges.csv' # input("Enter the file name containing the peak ranges (e.g., 'peak_ranges.csv'): ")
        peak_ranges = pd.read_csv(peak_range_file)  # Read the peak ranges from the specified file

        # print(peak_ranges)  # Display the peak ranges DataFrame for the user to see

    for file_name in data.keys():
        
        # print(f"Processing file: {file_name}")
        if manual_input == 'yes':
            amorphous_range = input("Enter the amorphous range (e.g., 20-30): ")
            amorphous_range = [float(x) for x in amorphous_range.split('-')]  # Convert the input range to a list of floats

            crystalline_range = input("Enter the crystalline range (e.g., 30-40): ")
            crystalline_range = [float(x) for x in crystalline_range.split('-')]  # Convert the input range to a list of floats

        else:
            omega = data[file_name]['omega']  # Get the omega value from the data dictionary
            file_index = peak_ranges.index[peak_ranges['Omega'] == omega]  # Find the index of the file in the peak ranges DataFrame

            amorphous_range = [float(x) for x in peak_ranges.loc[file_index, 'Amorphous Range'].values[0].split('-')]  # Get the amorphous range from the DataFrame
            crystalline_range = [float(x) for x in peak_ranges.loc[file_index, 'Crystalline Range'].values[0].split('-')]  # Get the crystalline range from the DataFrame

        # print("\n") # Print a new line for better readability

        data[file_name]['amorphous_range'] = amorphous_range
        data[file_name]['crystalline_range'] = crystalline_range
    
    return data

def Gaussian_function(current_x_value: float, amplitude, mean: float, standard_deviation: float, background: float): # Function to fit to
    """
    Calculates the Gaussian function value for a given x value.
    
    Parameters:
    current_x_value (float): The x value at which to evaluate the Gaussian function.
    
    amplitude (float): The amplitude of the Gaussian peak.
    
    mean (float): The mean (center) of the Gaussian peak.
    
    standard_deviation (float): The standard deviation of the Gaussian peak.

    background (float): The background value to be added to the Gaussian function.
    
    Returns:
    float: The value of the Gaussian function at the given x value.
    """

    func = amplitude*np.exp(-.5*((current_x_value-mean)/standard_deviation)**2)/(standard_deviation*np.sqrt(2*np.pi))+ background  # Calculate the Gaussian function value
    
    return func

def Peak_fit(data: dict): # Fit the peaks to a Gaussian function
    """
    Fits the peaks in the XRD data based on the provided ranges.
    
    Parameters:
    data (dict): A dictionary containing XRD data with file names as keys and 
    position/intensity as values, including amorphous and crystalline ranges.
    
    Returns:
    dict: An updated dictionary with fitted peak areas for each file.
    """

    for file_name, file_data in data.items(): # Iterate through each file in the data dictionary using key and value
        
        # print(f"Processing file: {file_name}")  # Print the file name being processed

        # Unpack the values from the file data
        position = file_data['Position']
        intensity = file_data['Intensity']
        amorphous_range_start_user, amorphous_range_end_user = file_data['amorphous_range']
        crystalline_range_start_user, crystalline_range_end_user = file_data['crystalline_range']

        # Step size calculation
        step_size = position[1] - position[0]  # Calculate the step size based on the position data

        # Convert the user input ranges to actual recorded values
        amorphous_range_start = position[(position >= amorphous_range_start_user) & (position <= amorphous_range_end_user)].min()
        amorphous_range_end = position[(position >= amorphous_range_start_user) & (position <= amorphous_range_end_user)].max()
        crystalline_range_start = position[(position >= crystalline_range_start_user) & (position <= crystalline_range_end_user)].min()
        crystalline_range_end = position[(position >= crystalline_range_start_user) & (position <= crystalline_range_end_user)].max()

        # Check range overlap, and set x values accordingly
        if (amorphous_range_start < crystalline_range_end and amorphous_range_end > crystalline_range_start):
            x_start = min(amorphous_range_start, crystalline_range_start) # Take the minimum of the two ranges
            x_end = max(amorphous_range_end, crystalline_range_end) # Take the maximum of the two ranges
            x_values = np.arange(x_start, x_end+0.5*step_size, step_size) # Create an array of x values for fitting

            x_values_amorphous = [i for i in x_values if i <=crystalline_range_start or i >= crystalline_range_end] # Get the x values for the amorphous range
            x_values_crystalline = [i for i in x_values if crystalline_range_start <= i <= crystalline_range_end] # Get the x values for the crystalline range
            

        else:
            x_values_amorphous = position[(position >= amorphous_range_start) & (position <= amorphous_range_end)]
            x_values_crystalline = position[(position >= crystalline_range_start) & (position <= crystalline_range_end)]
            
        # Determine y values for the Gaussian fit
        y_values_amorphous = intensity[(position >= amorphous_range_start) & (position <= amorphous_range_end) & ~((position >= crystalline_range_start) & (position <= crystalline_range_end))]  # Get the intensity values for the amorphous range
        y_values_crystalline = intensity[(position >= crystalline_range_start) & (position <= crystalline_range_end)]  # Get the intensity values for the crystalline range
        
        # Get guesses for the Gaussian parameters
        guess_background = intensity.min()  # Guess the background for the peaks

        guess_amplitude_amorphous = intensity[(position >= amorphous_range_start) & (position <= amorphous_range_end)].max() # Guess the amplitude for the amorphous peak
        guess_mean_amorphous = position[(position >= amorphous_range_start) & (position <= amorphous_range_end)].mean() # Guess the mean for the amorphous peak
        guess_std_amorphous = (amorphous_range_end - amorphous_range_start) / 6  # Guess the standard deviation for the amorphous peak
        guess_amorphous = [guess_amplitude_amorphous, guess_mean_amorphous, guess_std_amorphous, guess_background]  # Create a list of guesses for the amorphous peak

        guess_amplitude_crystalline = intensity[(position >= crystalline_range_start) & (position <= crystalline_range_end)].max() # Guess the amplitude for the crystalline peak
        guess_mean_crystalline = position[(position >= crystalline_range_start) & (position <= crystalline_range_end)].mean() # Guess the mean for the crystalline peak
        guess_std_crystalline = (crystalline_range_end - crystalline_range_start) / 6
        guess_crystalline = [guess_amplitude_crystalline, guess_mean_crystalline, guess_std_crystalline, guess_background]
        
        # Check size for fitting, and corrrect if necessary
        size_of_x_values_amorphous = len(x_values_amorphous)  # Get the size of the x values for the amorphous range
        size_of_y_values_amorphous = len(y_values_amorphous)  # Get the size of the y values for the amorphous range
        size_of_x_values_crystalline = len(x_values_crystalline)  # Get the size of the x values for the crystalline range
        size_of_y_values_crystalline = len(y_values_crystalline)  # Get the size of the y values for the crystalline range

        # Adjust the x and y values to ensure they are the same size for fitting
        if size_of_x_values_amorphous > size_of_y_values_amorphous:
            x_values_amorphous = x_values_amorphous[:size_of_y_values_amorphous] # Adjust the x values for the amorphous range to match the size of the y values
        elif size_of_x_values_amorphous < size_of_y_values_amorphous:
            y_values_amorphous = y_values_amorphous[:size_of_x_values_amorphous]

        if size_of_x_values_crystalline > size_of_y_values_crystalline:
            x_values_crystalline = x_values_crystalline[:size_of_y_values_crystalline]
        elif size_of_x_values_crystalline < size_of_y_values_crystalline:
            y_values_crystalline = y_values_crystalline[:size_of_x_values_crystalline]
        
        # Fit the Gaussian function to the data
        fit_amorphous, error_amorphous = curve_fit(Gaussian_function,x_values_amorphous,y_values_amorphous,p0=guess_amorphous)
        uncert_amorphous = np.sqrt(np.diag(error_amorphous))

        fit_crystalline, error_crystalline = curve_fit(Gaussian_function,x_values_crystalline,y_values_crystalline,p0=guess_crystalline)
        uncert_crystalline = np.sqrt(np.diag(error_crystalline))

        # Plot the fitted peaks for visualization
        plt.plot(x_values_amorphous, y_values_amorphous, 'ro', label='Amorphous Data')  # Plot the amorphous data points
        plt.plot(x_values_amorphous, Gaussian_function(x_values_amorphous, *fit_amorphous), 'b-', label='Amorphous Fit')  # Plot the fitted amorphous peak
        plt.plot(x_values_crystalline, y_values_crystalline, 'go', label='Crystalline Data')  # Plot the crystalline data points
        plt.plot(x_values_crystalline, Gaussian_function(x_values_crystalline, *fit_crystalline), 'r-', label='Crystalline Fit')  # Plot the fitted crystalline peak
        plt.xlabel('Position [°2θ]')  # Set the x-axis label
        plt.ylabel('Intensity [Counts]')  # Set the y-axis label
        plt.title(f'Peak Fitting for {file_name}')  # Set the title of the plot
        plt.legend()  # Show the legend
        plt.show()

        # Store the fitted parameters and uncertainties in the data dictionary
        data[file_name]['amorphous_fit'] = {
            'amplitude': fit_amorphous[0],
            'mean': fit_amorphous[1],
            'std_dev': fit_amorphous[2],
            'background': fit_amorphous[3],
            'uncertainty': uncert_amorphous,
            'amorphous_x_range_used': x_values_amorphous,  # Convert to list to save values for later use
            'amorphous_y_range_used': y_values_amorphous  # Convert to list to save values for later use
        }
        
        data[file_name]['crystalline_fit'] = {
            'amplitude': fit_crystalline[0],
            'mean': fit_crystalline[1],
            'std_dev': fit_crystalline[2],
            'background': fit_crystalline[3],
            'uncertainty': uncert_crystalline,
            'crystalline_x_range_used': x_values_crystalline,
            'crystalline_y_range_used': y_values_crystalline
        }
    
    return data

def Peak_Area_Calculation(data: dict):
    """
    Calculates the area under the fitted peaks for each file in the data dictionary.
    
    Parameters:
    data (dict): A dictionary containing XRD data with fitted peak parameters.
    
    Returns:
    dict: An updated dictionary with calculated peak areas for each file.
    """
    
    for file_name, file_data in data.items():
        # Calculate the area under the amorphous peak
        amplitude_amorphous = file_data['amorphous_fit']['amplitude']
        mean_amorphous = file_data['amorphous_fit']['mean']
        std_dev_amorphous = file_data['amorphous_fit']['std_dev']
        background_amorphous = file_data['amorphous_fit']['background']
        x_values_amorphous = file_data['amorphous_fit']['amorphous_x_range_used']
        area_amorphous = np.trapz(Gaussian_function(x_values_amorphous, amplitude_amorphous, mean_amorphous, std_dev_amorphous, background_amorphous), x=x_values_amorphous)  # Calculate the area under the amorphous peak

        # Calculate the area under the crystalline peak
        amplitude_crystalline = file_data['crystalline_fit']['amplitude']
        mean_crystalline = file_data['crystalline_fit']['mean']
        std_dev_crystalline = file_data['crystalline_fit']['std_dev']
        background_crystalline = file_data['crystalline_fit']['background']
        x_values_crystalline = file_data['crystalline_fit']['crystalline_x_range_used']
        area_crystalline = np.trapz(Gaussian_function(x_values_crystalline, amplitude_crystalline, mean_crystalline, std_dev_crystalline, background_crystalline), x=x_values_crystalline) # Calculate the area under the crystalline peak
        
        # Store the areas in the data dictionary
        data[file_name]['amorphous_area'] = area_amorphous
        data[file_name]['crystalline_area'] = area_crystalline
    
    return data

def Peak_Mixture_Calculation(data: dict): # Calculates the mixture fraction from amorphous and crystalline areas
    """
    Calculates the mixture of amorphous and crystalline areas for each file in the data dictionary.
    
    Parameters:
    data (dict): A dictionary containing XRD data with calculated peak areas.
    
    Returns:
    dict: An updated dictionary with calculated peak mixtures for each file.
    """
    
    for file_name, file_data in data.items():
        # Unpack needed values from data
        area_amorphous = file_data['amorphous_area']
        area_crystalline = file_data['crystalline_area']
        
        if area_amorphous > area_crystalline:
            crystalline_percentage =(area_crystalline/area_amorphous)
            amorphous_percentage = 1 - crystalline_percentage

        else:
            amorphous_percentage = (area_amorphous/area_crystalline)
            crystalline_percentage = 1 - amorphous_percentage
        
        data[file_name]['amorphous_percentage'] = amorphous_percentage
        data[file_name]['crystalline_percentage'] = crystalline_percentage
    
    return data

def Sample_Depth_Analysis(data: dict):
    """
    Analyzes the sample depth based on the calculated peak areas.
    
    Parameters:
    data (dict): A dictionary containing XRD data with calculated peak mixtures.
    
    Returns:
    dict: An updated dictionary with sample depth analysis results for each file.
    """
    
    def skin_depth(alpha: float): # Skin depth calculation (only used in this function)
        """
        Calculate the skin depth for a given attenuation coefficient.

        Parameters:
        alpha (float): Attenuation coefficient in m^-1.

        Returns:
        float: Skin depth in meters.
        """
        return 1 / alpha  # Skin depth in meters

    def real_depth(skin_depth: float, ω: float): # Real depth calculation (only used in this function)
        """
        Calculate the real depth for a given attenuation coefficient and angle.
        This calculates the depth at which 90% of the signal is attenuated.

        Parameters:
        skin_depth (float): Skin depth of mixture.
        ω (float): Angle in degrees.

        Returns:
        float: Real depth in micrometers.
        """
        return 5 * skin_depth * np.sin(np.radians(ω)) * 1e6  # Convert to micrometers
    
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

    for file_name, file_data in data.items():
        # Unpack needed values from data
        omega = file_data['omega']
        amorphous_percentage = file_data['amorphous_percentage']
        crystalline_percentage = file_data['crystalline_percentage']
        
        # website for index of refraction https://henke.lbl.gov/optical_constants/getdb2.html
        # Get indices of refraction based on the file name
        if 'p' in file_name:
            Pn = (1-6.96137067E-06)-(1.98567363E-07J)
        else:
            Bn = (1-6.94444225E-06)-(5.98799899E-09J)
        Sin = (1-7.57442876E-06)-(1.72761077E-07J)

        # Calculate the effective index of refraction based on the mixture
        effective_n = (amorphous_percentage * Pn + crystalline_percentage * Sin)

        # Calculate the attenuation coefficients
        effective_alpha = attenuation_coefficient(effective_n) # Attenuation coefficient in m^-1

        # Calculate skin depth and real depth
        skin = skin_depth(effective_alpha)
        effective_depth = real_depth(skin, omega)

        data[file_name]['effective_depth'] = effective_depth
    
    return data

def Store_Data(data: dict, path_name: str):
    """
    Stores the processed data into a CSV file.
    
    Parameters:
    data (dict): A dictionary containing processed XRD data.
    path_name (str): The path where the CSV file will be saved.
    """
    # print(data)  # Print the data dictionary to the console for debugging purposes
    
    df = pd.DataFrame.from_dict(data, orient='index')  # Convert the dictionary to a DataFrame
    df = df.sort_values(by='omega', ascending=True)  # Sort by Omega from smallest to largest
    df.to_csv(os.path.join(path_name, 'Processed_XRD_Data.csv'))  # Save the DataFrame to a CSV file

def bar_graph_depths():
    """
    Create a bar graph of the depths calculated for each mixture.
    """
    df = pd.read_csv('E://Full//p csv//Processed_XRD_Data.csv')
    x = np.arange(len(df['omega']))
    plt.bar(x, df['crystalline_percentage'] + df['amorphous_percentage'], label="Si")
    plt.bar(x, df['amorphous_percentage'], label="P")
    plt.xlabel('Omega (ω)')
    plt.ylabel('Mix')
    plt.title('Depths for Mixtures')
    plt.xticks(ticks=x, labels=[f"{df['omega'][i]}ω-{df['effective_depth'][i]:.2f}μ"for i in range(len(df['omega']))], rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.show()

def Peak_Area_Fitting(data_path: str):
    """
    Main function to perform peak area fitting on XRD data.
    
    Parameters:
    data_path (str): The path to the directory containing the XRD data files.
    """
    data_raw = XRD_Data_Dictionary(data_path)  # Create a dictionary with the XRD data
    data_ranges = Length_Peak(data_raw)  # Get the ranges for each file, only function that requires user input
    data_fit_parameters = Peak_fit(data_ranges)  # Fit the peaks to a Gaussian function
    data_peak_areas = Peak_Area_Calculation(data_fit_parameters)  # Calculate the area under the peaks
    data_mixture = Peak_Mixture_Calculation(data_peak_areas)  # Calculate the mixture of amorphous and crystalline areas
    data_depth = Sample_Depth_Analysis(data_mixture)  # Analyze the sample depth based on the calculated peak areas
    Store_Data(data_depth, data_path)  # Store the processed data into a CSV file

if __name__ == "__main__":

    # This assumes that for composite peaks, the crystalline peak is somewhere in the middle of the amorphous peak.

    Peak_Area_Fitting("E://Full//p csv")  # Call the main function with the path to the XRD data files

    bar_graph_depths() # Create a bar graph of the depths calculated for each mixture