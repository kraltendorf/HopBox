# **HopBox Inference**  #

Welcome to **HopBox Inference!** This is a simple, easy to use, Pytorch-Based FCN semantic segmentation model for to detect and extract Hop cone shape and color information from images. 

This model was trained on a dataset of about 12 images of hop cones, and is capable of detecting and segmenting Hop cones in images with a high degree of accuracy, then uses the segmentation to extract color and shape information. The model is too large to feature on GitHub so it needs to be downloaded from Google Drive using the publicly available Google Drive Folder [Publicly Available HopBox Resources](https://drive.google.com/drive/folders/137de5nLf2781ed__zGgdpxRoeVoNl9d4?usp=share_linkhttps://drive.google.com/drive/folders/137de5nLf2781ed__zGgdpxRoeVoNl9d4?usp=share_link)

For more information on the experimental setup, model architecture, color correction and results, please see the accompanying paper [paper].

This inference code is designed to be run on Google Colab or on any other platform with minimum adjustment. The code is designed to be run inference on a directory of images. The code will output three directories, one for **color corrected images**, another for **segmentation overlayed** on the original images, and a one more for the **segmented masks**. The code will also output a `.csv` file with the color and shape information for each Hop cone in the image.

To run the code, follow the steps below:

### **Step 1:** Prepare the working directory
In the working directory, create a folder that contains the images to run inference. The images can be in any format.
Make sure you have the following files in the working directory or in a subdirectory:
 - `HopBox_Infer.ipynb` # This is the inference code
 - `model.pth` # This is the model weights
 - `hopbox_functions.py` # Python file that contains all the functions used in the inference code
 - `ref_babel_avg_srgb.csv` # Contains the reference color values for the color correction (X-Rite ColorChecker Classic)

 If some files are in a subdirectory, make sure to change the path in the code to reflect the subdirectory.

 Open the `HopBox_Infer.ipynb` file in Google Colab to run inference.

### **Step 2:** Mount your Google Drive
Run the first cell to mount your Google Drive. This will allow you to access the files in your Google Drive from the Colab notebook.

```python
from google.colab import drive
drive.mount('/content/drive')
!ls "/content/drive/MyDrive/Colab notebook for Kayla" # Change the path to the directory where you have the files
```

This will also list the files in the directory you specified in the `!ls` command. Make sure the files are in the directory you specified.

### **Step 3:** Install dependencies
Run the second cell to install the dependencies. This will install the following packages:
 - `torch`      : Pytorch, used to load the model and run inference
 - `torchvision`: Pytorch vision, contains the model zoo
 - `torchaudio` : Pytorch audio, contains the audio models (optional)
 - `Pyzbar`     : Python library to read barcodes and QR codes
 - `libzbar0`   : C library to read barcodes and QR codes, required by Pyzbar for linux systems only

 ```python
 # install basic dependencies
!pip install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html
!apt install libzbar0 # if linux based machine
!pip install pyzbar

#Check current/allocated cuda GPU
!nvidia-smi
```

### **Step 4:** Load dependencies
Load the dependencies in the third cell. This will load useful packages and functions from the `hopbox_functions.py` file.

```python
from __future__ import print_function
from hopbox_functions import run_inference
import sys
import numpy as np
import os

sys.path.append(os.path.join(sys.path[0])) # add current directory to path
```

### **Step 5:** Load the Model and Data Paths
Load the model and the data paths in the fourth cell. Be sure to change the path to the model and the path to the directory containing the images to run inference on.

```python
Base_folder = "/content/drive/MyDrive/Colab notebook for Kayla/"
inference_dir = Base_folder + "NewImage_01/"

filess = os.listdir(inference_dir)
print(str(len(filess)) + " image files found")
print(filess)

Save_Images = 1 # 1 for saving images, and 0 for not saving images

model_dir = Base_folder + "best_models/hoptrial_250 epoch/hoptrial1.pth"

output_directory = Base_folder + "Results_DEC_22/NewImage_01"

if Save_Images == 1:
  mask_dir = output_directory + "/Masks/"
  label_dir = output_directory + "/Mask_labels/"
  corrected_dir = output_directory + "/Corrected_imgs/"

  try:
    os.makedirs(mask_dir)
    os.makedirs(label_dir)
    os.makedirs(corrected_dir)
  except:
    print("'" + output_directory+ "' already exists")
    
directories = [inference_dir,
               mask_dir,
               label_dir,
               corrected_dir,
               model_dir,
               output_directory
]
```
*`base_folder` is the path to the directory where you have the files. Change this to the path to the directory where you have the files.*

*This section also creates the output directories. If you do not want to save the images, set `Save_Images` to 0.*

*`inference_dir` is the path to the directory containing the images to run inference on. Change this to the path to the directory containing the images to run inference on.*

*`model_dir` is the path to the model weights. Change this to the path to the model weights.*

*`output_directory` is the path to the directory where you want to save the output. Change this to the path to the directory where you want to save the output.*

*`directories` is a list of all the directories. This is input to the `run_inference` function.*

### **Step 6:** Run inference all the images in the directory

Run the fifth cell to run inference on all the images in the directory. This will run inference on all the images in the directory and save the results in the output directories.

*`ref` is the reference color values for the color correction. This is a 24x3 array with the sRGB values for the reference colors. The reference color values are stored in the `ref_babel_avg_srgb.csv` file. If you want to use a different reference color chart, you can change the reference color values in the `ref_babel_avg_srgb.csv` file, or you can change the reference color values in the `load_csv_chart` function in the `hopbox_functions.py` file.*

*`scale` is the scale factor to resize the images. This is used to resize the images to the size of the network input. The default value is 1. If you want to resize the images, change the value of `scale` to the desired scale factor.*

```python
scale = 1; # should be based image size and network input size

# Reference Data 
def load_csv_chart(f):
    '''Input CSV's shape is (24, 3), with sRGB reference values
    '''
    data = np.loadtxt(f,delimiter=",")
    assert data.shape ==(24, 3)
    return data
ref = load_csv_chart(Base_folder + 'ref_babel_avg_srgb_d50.csv') # loads reference RGB color values

run_inference(ref=ref, directories=directories, filess=filess, scale=scale, Save_Images=Save_Images) # run inference on images in directory
```

_____
## **Samples of the input and output images are shown below:**

## Original Image vs Color Corrected Image

<img src="./ReadMe_files/original.png" />


## Hop Cone Model Detections

<img src= "./ReadMe_files/Detection.png" />

## Original vs color corrected Values

<img src="./ReadMe_files/Before_After_Color Correction.png" />

## .csv file output
<img src="./ReadMe_files/csv_.png" />


Citation:

```Journal Article (Peer Reviewed)```
```@article{10.1002/rob2.20202,```
```author = {John, Doe and Jane, Doe},```
```title = {HopBox: A Deep Learning Approach to Hop Cone Detection and Color Correction},```
```journal = {xxx},```
```volume = {xx},```
```number = {xx},```
```pages = {xx-xx},```
```year = {2020},```
```doi = {10.1002/rob2.20202},```
```url = {https://doi.org/10.1002/rob2.20202},```
```eprint = {https://doi.org/10.1002/rob2.20202}}```
