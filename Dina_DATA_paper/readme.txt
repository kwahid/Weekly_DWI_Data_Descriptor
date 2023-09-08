run fitting_script for ADC map generation.

Beware of the complication of DICOM format from different MR scanners, please use https://github.com/rordenlab/dcm2niix to convert DICOM to nifti format.


data_order options: {'xynz','xyzn','xyzfile'}, for 4D DWI file in nifti format choose 'xyzn'