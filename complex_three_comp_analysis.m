function complex_three_comp_analysis(mag_data_file_name, phase_data_file_name, cls_file_name)

% *************************************************************************
% function complex_three_comp_analysis(mag_data_file_name, phase_data_file_name, cls_file_name)
%
% AUTHOR: Eva Alonso Ortiz, PhD 
% email: eva.alonso.ortiz@gmail.com
%
% DESCRIPTION: Complex 3-component analysis of multi-echo T2* decay. Please refer to 
% NeuroImage 182 (2018) 370?378 for details.
%
% DATE: Feb 2019
% 
%**************************************************************************

% =========================== Header ==================================== %
this_fname = 'complex_three_comp_analysis';
this_info = sprintf('%-20s : ',this_fname);
fprintf([this_info, 'Current date and time: %s\n'], datestr(now));
% =========================================================================

B0 = 3;

% Fitting average ROI data? Yes=1/No=0
ROI_flag = 1;

%%-------------------------------------------------------------------------
%% check existence of data files
%%-------------------------------------------------------------------------

[fid, message] = fopen(mag_data_file_name,'r');
if(fid == -1)
  error(sprintf('\nError in multi_comp_fit: cannot find input data file %s\n', mag_data_file_name));
else
  fclose(fid);
end

[fid, message] = fopen(phase_data_file_name,'r');
if(fid == -1)
  error(sprintf('\nError in multi_comp_fit: cannot find input data file %s\n', phase_data_file_name));
else
  fclose(fid);
end


% check existence of mask file
[fid, message] = fopen(cls_file_name,'r');
if(fid == -1)
    error(sprintf('\nError in multi_comp_fit: cannot find input mask file %s\n', cls_file_name));
else
    fclose(fid);
    val = input('Enter the tissue flag to be used for processing (1->CSF, 2->GM, 3->WM): ');
    switch val
        case 1
            tissue_flag = 1;
        case 2
            tissue_flag = 2;
        case 3 
            tissue_flag = 3;
        otherwise
            warning(sprintf('Invalid tissue flag entered. Proceeding analysis with no mask.'))
    end                        
end



%%-------------------------------------------------------------------------
%% open data file
%%-------------------------------------------------------------------------

[mag_data_desc,mag_data_vol] = niak_read_minc(mag_data_file_name);
[phase_data_desc,phase_data_vol] = niak_read_minc(phase_data_file_name);

mag_data_dim = mag_data_desc.info.dimensions;
mag_data_slices = mag_data_dim(1,3);
mag_data_height = mag_data_dim(1,1);
mag_data_width = mag_data_dim(1,2);
mag_data_voxels = mag_data_height*mag_data_width;
num_echoes = mag_data_dim(1,4);

%%-------------------------------------------------------------------------
%% open mask files
%%-------------------------------------------------------------------------

% open tissue classification mask
[mask_desc,mask_vol] = niak_read_minc(cls_file_name);

mask_dim = mask_desc.info.dimensions;
mask_voxels = mask_dim(1,1)*mask_dim(1,2);

% check that mask and data_vol are the same dimensions
if mask_voxels ~= mag_data_voxels
    error(sprintf('\nError in multi_comp_fit: Tissue mask file dimensions do not match data image file.\n')); 
end

%%-------------------------------------------------------------------------
%% caculate echo times
%%-------------------------------------------------------------------------

echo_times = calc_echo_times(num_echoes);
echo_times = 1e-3*echo_times;

mag_data_vol = mag_data_vol(:,:,:,1:num_echoes);
phase_data_vol = phase_data_vol(:,:,:,1:num_echoes);

%%-------------------------------------------------------------------------
%% caculate real and imaginary volumes
%%-------------------------------------------------------------------------

real_data_vol = zeros(mag_data_height,mag_data_width,mag_data_slices,num_echoes);
imag_data_vol = zeros(mag_data_height,mag_data_width,mag_data_slices,num_echoes);

for slice = 1:mag_data_slices

    for i = 1:mag_data_height
        for j = 1:mag_data_width

            real_data_vol(i,j,slice,:) = mag_data_vol(i,j,slice,:).*cos(phase_data_vol(i,j,slice,:));
            imag_data_vol(i,j,slice,:) = mag_data_vol(i,j,slice,:).*sin(phase_data_vol(i,j,slice,:));

        end
    end
end

real_data_vol(isnan(real_data_vol)) = 0 ;
imag_data_vol(isnan(imag_data_vol)) = 0 ;

%%-------------------------------------------------------------------------
%% data fitting and analysis
%%------------------------------------------------------------------------- 

if (ROI_flag == 1)
    
    for echo = 1:num_echoes  
        mag_roi_data = squeeze(mag_data_vol(:,:,1,echo)).*mask_vol(:,:);    
        mag_roi_data = (reshape(mag_roi_data, mag_data_voxels, 1));
        mean_mag_roi_data(echo) = mean(nonzeros(mag_roi_data));
        
        real_roi_data = squeeze(real_data_vol(:,:,1,echo)).*mask_vol(:,:);    
        real_roi_data = (reshape(real_roi_data, mag_data_voxels, 1));
        mean_real_roi_data(echo) = mean(nonzeros(real_roi_data));
        
        imag_roi_data = squeeze(imag_data_vol(:,:,1,echo)).*mask_vol(:,:);    
        imag_roi_data = (reshape(imag_roi_data, mag_data_voxels, 1));
        mean_imag_roi_data(echo) = mean(nonzeros(imag_roi_data));
    end
    
    %-----------------------------------------------
    % Do 3-component multi-exponential fitting
    %-----------------------------------------------  
    [RE_fit_multi, IM_fit_multi] = complex_3_comp_fit(mean_real_roi_data, mean_imag_roi_data, echo_times, B0);
    [A_fit, T2s_fit, resnorm] = three_comp_fit(mean_mag_roi_data, echo_times,B0);

    % Calculate complex 3CF parameters
    complex_3CF_A_MW = (RE_fit_multi(1)+IM_fit_multi(1))/2;
    complex_3CF_T2s_MW = (RE_fit_multi(2)+IM_fit_multi(2))/2; 
    complex_3CF_omega_MW = (RE_fit_multi(3)+IM_fit_multi(3))/2;

    complex_3CF_A_MW_RE = RE_fit_multi(1);
    complex_3CF_T2s_MW_RE = RE_fit_multi(2); 
    complex_3CF_omega_MW_RE = RE_fit_multi(3);

    complex_3CF_delf_MW = complex_3CF_omega_MW/(2*pi);
    complex_3CF_delf_MW_RE = complex_3CF_omega_MW_RE/(2*pi);

    complex_3CF_A_EW = RE_fit_multi(4);
    complex_3CF_T2s_EW = RE_fit_multi(5); 

    complex_3CF_A_AW = (RE_fit_multi(6)+IM_fit_multi(4))/2;
    complex_3CF_T2s_AW = (RE_fit_multi(7)+IM_fit_multi(5))/2;
    complex_3CF_omega_AW = (RE_fit_multi(8)+IM_fit_multi(6))/2;

    complex_3CF_A_AW_RE = RE_fit_multi(6);
    complex_3CF_T2s_AW_RE = RE_fit_multi(7);
    complex_3CF_omega_AW_RE = RE_fit_multi(8);

    complex_3CF_delf_AW = complex_3CF_omega_AW/(2*pi);
    complex_3CF_delf_AW_RE = complex_3CF_omega_AW_RE/(2*pi);

    MWF_complex_3CF = complex_3CF_A_MW/(complex_3CF_A_MW+complex_3CF_A_AW+complex_3CF_A_EW) ;
    MWF_complex_3CF_RE = complex_3CF_A_MW_RE/(complex_3CF_A_MW_RE+complex_3CF_A_AW_RE+complex_3CF_A_EW) ;
    
    AWF_complex_3CF = complex_3CF_A_AW/(complex_3CF_A_MW+complex_3CF_A_AW+complex_3CF_A_EW) ;
    AWF_complex_3CF_RE = complex_3CF_A_AW_RE/(complex_3CF_A_MW_RE+complex_3CF_A_AW_RE+complex_3CF_A_EW) ;
    
    EWF_complex_3CF = complex_3CF_A_EW/(complex_3CF_A_MW+complex_3CF_A_AW+complex_3CF_A_EW) ;
    EWF_complex_3CF_RE = complex_3CF_A_EW/(complex_3CF_A_MW_RE+complex_3CF_A_AW_RE+complex_3CF_A_EW) ;
    
    MWF_3CMF = A_fit.MW/(A_fit.MW+A_fit.AW+A_fit.EW);
    AWF_3CMF = A_fit.AW/(A_fit.MW+A_fit.AW+A_fit.EW);
    EWF_3CMF = A_fit.EW/(A_fit.MW+A_fit.AW+A_fit.EW);
    
    threeCMF_T2s_MW = T2s_fit.MW;
    threeCMF_T2s_AW = T2s_fit.AW;
    threeCMF_T2s_EW = T2s_fit.EW;

          
    mean_roi_data = [MWF_complex_3CF_RE AWF_complex_3CF_RE EWF_complex_3CF_RE complex_3CF_T2s_MW complex_3CF_T2s_AW complex_3CF_T2s_EW complex_3CF_delf_MW_RE complex_3CF_delf_AW_RE MWF_3CMF AWF_3CMF EWF_3CMF threeCMF_T2s_MW threeCMF_T2s_AW threeCMF_T2s_EW];
    dlmwrite('3CCF.dat',[mean_roi_data], '-append');

else

    % initialize all data vecotrs to zero
    complex_3CF_A_MW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_MW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_omega_MW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_A_EW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_EW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_A_AW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_AW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_omega_AW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_delf_MW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_delf_AW = zeros(mag_data_height,mag_data_width,mag_data_slices);
    MWF_complex_3CF = zeros(mag_data_height,mag_data_width,mag_data_slices);
    EWF_complex_3CF = zeros(mag_data_height,mag_data_width,mag_data_slices);
    AWF_complex_3CF = zeros(mag_data_height,mag_data_width,mag_data_slices);
    
    complex_3CF_A_MW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_MW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_omega_MW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_A_EW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_EW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_A_AW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_T2s_AW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_omega_AW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_delf_MW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    complex_3CF_delf_AW_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    MWF_complex_3CF_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    EWF_complex_3CF_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);
    AWF_complex_3CF_RE = zeros(mag_data_height,mag_data_width,mag_data_slices);

    for slice = 1:mag_data_slices

        for i = 1:mag_data_height
            for j = 1:mag_data_width

                % only fill in what corresponds to the mask
                if mask_vol(i,j,slice) == tissue_flag

                    %-----------------------------------------------
                    % Do 3-component multi-exponential fitting
                    %-----------------------------------------------  
 
                    [RE_fit_multi, IM_fit_multi] = complex_3_comp_fit(squeeze(real_data_vol(i,j,slice,:))', squeeze(imag_data_vol(i,j,slice,:))', echo_times, B0);
                    
                    % Calculate complex 3CF parameters
                    complex_3CF_A_MW(i,j,slice) = (RE_fit_multi(1)+IM_fit_multi(1))/2;
                    complex_3CF_T2s_MW(i,j,slice) = (RE_fit_multi(2)+IM_fit_multi(2))/2; 
                    complex_3CF_omega_MW(i,j,slice) = (RE_fit_multi(3)+IM_fit_multi(3))/2;

                    complex_3CF_A_MW_RE(i,j,slice) = RE_fit_multi(1);
                    complex_3CF_T2s_MW_RE(i,j,slice) = RE_fit_multi(2); 
                    complex_3CF_omega_MW_RE(i,j,slice) = RE_fit_multi(3);

                    complex_3CF_delf_MW(i,j,slice) = complex_3CF_omega_MW(i,j,slice)/(2*pi);
                    complex_3CF_delf_MW_RE(i,j,slice) = complex_3CF_omega_MW_RE(i,j,slice)/(2*pi);

                    complex_3CF_A_EW_RE(i,j,slice) = RE_fit_multi(4);
                    complex_3CF_T2s_EW_RE(i,j,slice) = RE_fit_multi(5); 

                    complex_3CF_A_AW(i,j,slice) = (RE_fit_multi(6)+IM_fit_multi(4))/2;
                    complex_3CF_T2s_AW(i,j,slice) = (RE_fit_multi(7)+IM_fit_multi(5))/2;
                    complex_3CF_omega_AW(i,j,slice) = (RE_fit_multi(8)+IM_fit_multi(6))/2;

                    complex_3CF_A_AW_RE(i,j,slice) = RE_fit_multi(6);
                    complex_3CF_T2s_AW_RE(i,j,slice) = RE_fit_multi(7);
                    complex_3CF_omega_AW_RE(i,j,slice) = RE_fit_multi(8);

                    complex_3CF_delf_AW(i,j,slice) = complex_3CF_omega_AW(i,j,slice)/(2*pi);
                    complex_3CF_delf_AW_RE(i,j,slice) = complex_3CF_omega_AW_RE(i,j,slice)/(2*pi);

                    MWF_complex_3CF(i,j,slice) = complex_3CF_A_MW(i,j,slice)/(complex_3CF_A_MW(i,j,slice)+complex_3CF_A_AW(i,j,slice)+complex_3CF_A_EW(i,j,slice)) ;
                    AWF_complex_3CF(i,j,slice) = complex_3CF_A_AW(i,j,slice)/(complex_3CF_A_MW(i,j,slice)+complex_3CF_A_AW(i,j,slice)+complex_3CF_A_EW(i,j,slice)) ;
                    EWF_complex_3CF(i,j,slice) = complex_3CF_A_EW(i,j,slice)/(complex_3CF_A_MW(i,j,slice)+complex_3CF_A_AW(i,j,slice)+complex_3CF_A_EW(i,j,slice)) ;

                    MWF_complex_3CF_RE(i,j,slice) = complex_3CF_A_MW_RE(i,j,slice)/(complex_3CF_A_MW_RE(i,j,slice)+complex_3CF_A_AW_RE(i,j,slice)+complex_3CF_A_EW_RE(i,j,slice)) ;
                    AWF_complex_3CF_RE(i,j,slice) = complex_3CF_A_AW_RE(i,j,slice)/(complex_3CF_A_MW_RE(i,j,slice)+complex_3CF_A_AW_RE(i,j,slice)+complex_3CF_A_EW_RE(i,j,slice)) ;
                    EWF_complex_3CF_RE(i,j,slice) = complex_3CF_A_EW_RE(i,j,slice)/(complex_3CF_A_MW_RE(i,j,slice)+complex_3CF_A_AW_RE(i,j,slice)+complex_3CF_A_EW_RE(i,j,slice)) ;

                end

            end
        end

    end

    % Create  maps 
    mask_desc.info.dimensions(1,4) = 1;
    file_name = 'mwf_3CCF.mnc';
    data = 100*MWF_complex_3CF;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'awf_3CCF.mnc';
    data = 100*AWF_complex_3CF;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'ewf_3CCF.mnc';
    data = 100*EWF_complex_3CF;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'delf_AW_3CCF.mnc';
    data = complex_3CF_delf_AW;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'delf_MW_3CCF.mnc';
    data = complex_3CF_delf_MW;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    mask_desc.info.dimensions(1,4) = 1;
    file_name = 'mwf_3CCF_RE.mnc';
    data = 100*MWF_complex_3CF_RE;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'awf_3CCF_RE.mnc';
    data = 100*AWF_complex_3CF_RE;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'ewf_3CCF_RE.mnc';
    data = 100*EWF_complex_3CF_RE;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'delf_AW_3CCF_RE.mnc';
    data = complex_3CF_delf_AW_RE;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

    file_name = 'delf_MW_3CCF_RE.mnc';
    data = complex_3CF_delf_MW_RE;
    mask_desc.file_name = file_name;
    niak_write_minc(mask_desc,data);

end


 end


