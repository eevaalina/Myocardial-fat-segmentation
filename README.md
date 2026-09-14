# Myocardial-fat-segmentation
Python pipeline for automated detection of fat particles in the cardiac myocardium from CT images.

The pipeline processes DICOM CT data, converts the selected CT series to NIfTI, performs cardiac segmentation using TotalSegmentator, detects fat-density particles within the myocardium, and exports fat particle measurements as Excel files.
This code has been used in the study: URL coming…

# OVERVIEW
The analysis pipeline consists of the following steps:
1.	Search patient folders for DICOM CT series
2.	Select the CT series with the largest number of DICOM images
3.	Convert the selected DICOM series to NIfTI
4.	Perform cardiac segmentation using the TotalSegmentator "heartchambers_highres" model
5.	Fat particles are identified within the myocardial mask using a CT attenuation range of −200 to −40 HU
6.	Calculate particle volume and mean HU
7.	Detect particles with possible metal artifacts 
8.	Calculate noise-control measurements from the aorta area
9.	Save visualization images and Excel files for fat-particle and noise-control assessment

# Requirements

This code has been written for Python 3.10. and the code uses the following Python packages:

•	NumPy <br>
•	SimpleITK<br>
•	NiBabel<br>
•	SciPy<br>
•	Matplotlib<br>
•	pandas<br>
•	openpyxl<br>
Install the required packages: <br>
`pip install numpy SimpleITK nibabel scipy matplotlib pandas openpyxl`


# Configuration
Before running the analysis, configure the input and output paths at the beginning of the Python script.
The following paths must be configured:
`dicom_root = r"/path/to/DICOM"`
dicom_root should point to the main directory containing the patient folders and images are expected to be as ‘.dcm’ .

`noise_folder = r"/path/to/noise_results"`
noise_folder specifies where the aortic noise-control Excel files are saved.

`results_folder = r"/path/to/results"`
results_folder specifies where the particle-analysis Excel files are saved.

The analysis can be started from a specific patient index by changing:
`start_index = 0`
For example, setting start_index = 10 starts the analysis from patient index 10.

Input data:

The pipeline expects patient folders containing DICOM CT data.
Example structure:
<img width="524" height="264" alt="image" src="https://github.com/user-attachments/assets/a7fe4090-4d19-4827-8a7f-1442f7db4bff" />

CT series selection:

The script searches for DICOM series within each patient folder and the series containing the largest number of DICOM images is selected for further analysis.
The selected series is then converted from DICOM to NIfTI using SimpleITK.

# Segmentation
The pipeline uses TotalSegmentator for automated cardiac segmentation with the heartchambers_highres task.
Install and configure TotalSegmentator according to the official documentation:
https://github.com/wasserth/TotalSegmentator

Cardiac segmentation:

TotalSegmentator is used with:
`task = "heartchambers_highres"` <br>
The following anatomical masks are retained: <br>
•	Aorta <br>
•	Myocardium <br>
•	Left ventricle <br>
These masks are used for subsequent fat-particle analysis and visualization.

Fat particle detection:

Fat-density candidates are identified within the myocardial mask using the following CT attenuation range: -200 HU ≤ HU ≤ -40 HU
This identifies regions within the myocardium with fat-like CT attenuation.

Edge-particle filtering:

To reduce false-positive particles originating from the boundary of the myocardium and left ventricle, the myocardium and left-ventricle masks are combined, and the combined mask is then eroded using a spherical structuring element with radius 3 voxels. Particles that touch the resulting outer-edge region are removed as possible epicardial fat.

Individual particles are identified using connected-component labeling.
For every retained particle, the following measurements are calculated:

<img width="531" height="157" alt="image" src="https://github.com/user-attachments/assets/e9752f82-9b42-4647-a8b1-385971935991" />


Metal-artifact detection: 

A metal mask is generated using: HU ≥ 3000
A particle receives a positive metal flag if a voxel with HU ≥ 3000 is present on any CT slice containing that particle. The metal flag is used for quality control, as metal artifacts may cause false-positive fat detections, and all image stacks with a metal flag are recommended to be visually evaluated.

Noise control :

Noise control is performed using the aorta.
The aorta mask is first eroded using the same spherical structuring element with radius 3 to reduce the contribution of boundary voxels.
The CT values within the resulting aortic mask are used to calculate:
•	Minimum HU
•	Maximum HU
•	Mean HU
•	Standard deviation
•	SNR
The noise-control results are saved separately as an Excel file.

# Output
For each CT slice, the pipeline saves a two-panel visualization as a .png to visually assess the quality of segmentations.
1. Original CT with aorta, myocardium, left ventricle, and outer-edge masks overlaid.
2. Original CT with detected fat particles overlaid.

For each patient, the converted CT image and the selected cardiac segmentation masks are stored in a MASK folder inside the corresponding patient directory. The converted CT image is saved as a NIfTI file named according to the patient index. 
The segmentation results are stored in a separate mask folder and include the aorta, myocardium, and left ventricle masks. The pipeline also saves a two-panel PNG visualization for each CT slice to enable the assessment of segmentations. 
The pipeline saves the noise-control measurements as a separate Excel file in the configured noise output folder. The particle-level analysis results, including particle volume, number of voxels, mean HU, and metal-artifact flag, are saved as a separate Excel file in the configured results output folder. 
The exact filenames depend on the patient index.

Results: 

The particle results are exported to Excel.
The result file contains CT acquisition metadata followed by the particle-level measurements. The Metal Flag is set to True when possible metal artefacts are detected on the same CT slice as the particle. All patients with possible metal artefacts should be visually reviewed by a qualified professional.

Example from results file structure:

<img width="604" height="227" alt="image" src="https://github.com/user-attachments/assets/6edfcdfd-bc77-4e8d-91c5-71702f739b47" />


# License

This project depends on external software, which has its own licensing terms. Please check the official TotalSegmentator website for its current license and usage conditions: https://github.com/wasserth/TotalSegmentator

# Acknowledgements

This project uses TotalSegmentator for automated anatomical segmentation.
TotalSegmentator: https://github.com/wasserth/TotalSegmentator




