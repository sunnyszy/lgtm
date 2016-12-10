%%
% The MIT License (MIT)
% Copyright (c) 2016 Ethan Gaebel <egaebel@vt.edu>
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.
%

%% Runs the SpotFi test over the passed in data files which each contain CSI data for many packets
function output_top_aoa = spotfi(sanitized_csi)
    % Acquire smoothed CSI matrix
    smoothed_sanitized_csi = smooth_csi(sanitized_csi);
    % Run SpotFi's AoA-ToF MUSIC algorithm on the smoothed and sanitized CSI matrix
    output_top_aoa = aoa_tof_music(smoothed_sanitized_csi);
end


function [estimated_aoa] = aoa_tof_music(x)
    %% DEBUG AND OUTPUT VARIABLES-----------------------------------------------------------------%%
    % Flow Variables 
    global OUTDOOR_FLAG
    
    % Output Variables
    global OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH

    % Phisical constant
    global n_subcarrier;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Data covarivance matrix
    R = x * x'; 
    % Find the eigenvalues and eigenvectors of the covariance matrix
    [eigenvectors, eigenvalue_matrix] = eig(R);
    
    if ~issorted(diag(eigenvalue_matrix))
        [eigenvalue_matrix,I] = sort(diag(eigenvalue_matrix));
        eigenvectors = eigenvectors(:, I);
    end
    
    % Find max eigenvalue for normalization
    max_eigenvalue = eigenvalue_matrix(n_subcarrier,n_subcarrier);

    for ii = 1:size(eigenvalue_matrix, 1)
        eigenvalue_matrix(ii, ii) = eigenvalue_matrix(ii, ii) / max_eigenvalue;
    end
    
    % Find the largest decrease ratio that occurs between the last 10 elements (largest 10 elements)
    % and is not the first decrease (from the largest eigenvalue to the next largest)
    % Compute the decrease factors between each adjacent pair of elements, except the first decrease
    start_index = size(eigenvalue_matrix, 1) - 2;
    end_index = start_index - 10;
    decrease_ratios = zeros(start_index - end_index + 1, 1);
    k = 1;
    for ii = start_index:-1:end_index
        temp_decrease_ratio = eigenvalue_matrix(ii + 1, ii + 1) / eigenvalue_matrix(ii, ii);
        decrease_ratios(k, 1) = temp_decrease_ratio;
        k = k + 1;
    end
    [~, max_decrease_ratio_index] = max(decrease_ratios);

    index_in_eigenvalues = size(eigenvalue_matrix, 1) - max_decrease_ratio_index;
    num_computed_paths = size(eigenvalue_matrix, 1) - index_in_eigenvalues + 1;
    
    % Estimate noise subspace
    column_indices = 1:(size(eigenvalue_matrix, 1) - num_computed_paths);
    eigenvectors = eigenvectors(:, column_indices); 
    % Angle in degrees 
    theta = 0:1:180; 
    % time in seconds
    %
    % TODO: put outdoor variable
    if ~OUTDOOR_FLAG
        %tau = 0:(0.2 * 10^-9):(20 * 10^-9); % 0 - 6m
        tau = 0:(1.0 * 10^-9):(200 * 10^-9);
    else
        tau = 0:(2 * 10^-9):(200 * 10^-9); % 0 - 60m
    end
    Pmusic = zeros(length(theta), length(tau));
    % Angle of Arrival Loop (AoA)
    for ii = 1:length(theta)
        % Time of Flight Loop (ToF)
        for jj = 1:length(tau)
            steering_vector = compute_steering_vector(theta(ii), tau(jj));
            PP = steering_vector' * (eigenvectors * eigenvectors') * steering_vector;
            Pmusic(ii, jj) = abs(1 /  PP);
        end
    end

    if OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH 
        % Theta (AoA) & Tau (ToF) 3D Plot
        figure('Name', 'AoA & ToF MUSIC Peaks', 'NumberTitle', 'off')
        mesh(tau, theta, Pmusic)
        xlabel('Time of Flight')
        ylabel('Angle of Arrival in degrees')
        zlabel('Spectrum Peaks')
        title('AoA and ToF Estimation from Modified MUSIC Algorithm')
        grid on
    end
    
    [~, ind] = max(Pmusic(:));
    [aoa_ind, ~] = ind2sub([length(theta), length(tau)], ind);
    estimated_aoa = theta(aoa_ind);


end

%% Computes the steering vector for SpotFi. 
function steering_vector = compute_steering_vector(theta, tau)

global delta_f channel_frequency d c subset_x subset_y;
    
    n_row = subset_x*subset_y;
    
    help_tau = repmat(1:subset_x, [subset_y, 1]) - 1;
    help_theta = repmat(transpose(1:subset_y), [1, subset_x]) - 1;
    
    steering_matrix = (omega_tof_phase(tau).^help_tau).*(phi_aoa_phase(theta).^help_theta);
    steering_vector = reshape(transpose(steering_matrix), [n_row, 1]);
    
function omega = omega_tof_phase(tau) 
    omega = exp(-1i * 2 * pi * delta_f * tau);
end

function phi = phi_aoa_phase(theta)
    phi = exp(-1i * 2 * pi * d * cos(deg2rad(theta)) * (channel_frequency / c));
end

end


function smoothed_csi = smooth_csi(csi)
    global n_subcarrier n_antenna subset_x subset_y
    n_row = subset_x*subset_y;
    n_column = (n_subcarrier-subset_x+1)*(n_antenna-subset_y+1);
    smoothed_csi = zeros(n_row, n_column);
    for i = 1:n_column
        % index of subcarrier and antenna
        [ind_s, ind_a] = ind2sub([n_subcarrier-subset_x+1, n_antenna-subset_y+1], i);
%         disp([ind_a, ind_a+subset_y-1, ind_s, ind_s+subset_x-1]);
        tmp_matrix = csi(ind_a:ind_a+subset_y-1, ind_s:ind_s+subset_x-1);
        smoothed_csi(:,i) = reshape(transpose(tmp_matrix), [n_row, 1]);
    end  
end