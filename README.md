# HopBox
 The HopBox is an a simply constructed lightbox and camera mount with an associated image analysis pipeline for evaluating hop cone morphology in the context of hop breeding or research programs. To create the HopBox we used readily available materials from the local hardware store (see Resources). To develop the image analysis pipeline, we used a fully convolutional network to train a model using 12 manually masked images (see Additonal Online Resources). The pipeline was evaluated by imaging 500 hop cones each from 15 experimental hop genotypes. Using the output from the HopBox, we were able to calculate genotype means, identify significant differences between genotypes for selection purposes, calculate heritability, correlations between traits, and conducted several analyses to determine optimal sample sizes of cones (see Output and Scripts). The pipeline can be conveniently run using the Google Colab platform (see Scripts).
 
This work has been submitted for publication in the Plant Phenome Journal.

## This respository contains: 

### Data
image weights, dry matter data

### Output
all output produced by the pipeline

### Resources
materials list and instructions for making your own hop box using readily available materials

### Scripts
google colab code for cone detection (please read the README), color correction, barcode reading, and data output pipeline, as well as various r scripts for analyzing the data, running analyses of variance, calculating heritability, creating figures, and optimizing cone sample numbers.

### Additional Online Resources
please visit our [Publicly Available HopBox Resources Folder](https://drive.google.com/drive/folders/137de5nLf2781ed__zGgdpxRoeVoNl9d4?usp=share_link) to download the model (too big for GitHub), access the training and test images that were used, and to download a set of images to test the pipeline. 
