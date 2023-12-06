import os
import re
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from tqdm import tqdm  # Import tqdm for progress bar

# Subdirectory path containing data files
directory_path = "data"  # Replace with the actual subdirectory path

# Get a list of all files in the directory
file_list = sorted([f for f in os.listdir(directory_path) if f.endswith('.txt')], key=lambda x: int(re.search(r'\d+', x).group()))

# Create a list to store images
images = []

# Iterate through each file using tqdm for a progress bar
for file_name in tqdm(file_list, desc="Processing Files", unit="file"):
    # Construct the full file path
    file_path = os.path.join(directory_path, file_name)

    # Read data from the file
    data_array = np.loadtxt(file_path)

    # Get the shape of the array
    rows, cols = data_array.shape

    # Plot grayscale image
    plt.imshow(data_array, cmap='gray_r', vmin=0, vmax=1, aspect='equal', extent=[0, cols, 0, rows])
    plt.colorbar()

    # Extract the number from the file name
    file_number = int(re.search(r'\d+', file_name).group())

    # Set title and labels
    plt.title(f'Grayscale Image - {file_number}')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    # Set axis ticks with an interval of 2
    plt.xticks(np.arange(0, cols+1, 5))
    plt.yticks(np.arange(0, rows+1, 5))

    # Save the image and add it to the list
    output_file = f"{file_number}_plot.png"
    plt.savefig(output_file)
    images.append(Image.open(output_file))

    # Close the current figure
    plt.close()

    # Delete intermediate image file
    os.remove(output_file)

# Save images as a GIF animation with a duration of 100 milliseconds
gif_output_file = "./plot/result.gif"
images[0].save(gif_output_file, save_all=True, append_images=images[1:], duration=100, loop=0)

# Save the last image separately
last_image_file = "./plot/final.png"
images[-1].save(last_image_file)

# Display the last image
images[-1].show()