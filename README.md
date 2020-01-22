# MWI
Myelin Water Imaging

- complex_three_comp_analysis.m

    DESCRIPTION:

    complex_three_comp_analysis.m is used with magnitude and phase multi-gradient echo images and a tissue mask (to limit     
    processing times) to generate myelin water fraction (MWF) maps based on 3-component complex fitting and 3-component 
    magnitude fitting (see NeuroImage 182 (2018) 370-378 for details).

    USAGE:

    complex_three_comp_analysis(mag_data_file_name, phase_data_file_name, cls_file_name)

    mag_data_file_name : multi-gradient echo magnitude image in minc format
    phase_data_file_name : multi-gradient echo phase image in minc format
    cls_file_name : tissue mask to be used for processing (1->CSF, 2->GM, 3->WM) in minc format

    NOTES:

    - Tested under Matlab 2015b
    - niak (https://www.nitrc.org/projects/niak/) must be downloaded and added to the matlab path

    CITATION

    If you use this code in your work, please cite:
    
    Alonso-Ortiz, E., Levesque, I. R., and Pike, G. B. Impact of magnetic susceptibility anisotropy at 3 T and 7 T on T2*-
    based myelin water fraction imaging. Neuroimage 182:370-378 (2018).

    
