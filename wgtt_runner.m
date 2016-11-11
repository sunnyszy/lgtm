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

function [top_aoas] = wgtt_runner()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    globals_init()
    
    %flow control
    global SIMULATION SIMULAIION_ALWAYS_GENERATE_DATA AOA_EST_MODE
    SIMULATION = false;
    SIMULAIION_ALWAYS_GENERATE_DATA = true;
    AOA_EST_MODE = 'MUSIC';%'MUSIC' 'SPOTFI';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the full path to the currently executing file and change the
    % pwd to the folder this file is contained in...
    [current_directory, ~, ~] = fileparts(mfilename('fullpath'));
    cd(current_directory);
    % Paths for the csitool functions provided
    path('./Atheros_csi', path);
    
    %% main
    if SIMULATION
        name_base = 'simulation_tmp';
        if ~exist('test-data/simulation_tmp.mat') || SIMULAIION_ALWAYS_GENERATE_DATA
            generate_simulation_data(['test-data/' name_base], 4000);
        end
    else   
        name_base = '135';
    end

    data_file = ['test-data/' name_base];
    top_aoas = run(data_file);
    
    save(['test-output/' name_base]);
end


%% Runs the SpotFi test over the passed in data files which each contain CSI data for many packets
% data_files -- a cell array of file paths to data files
function output_top_aoa = run(data_file)
    %% DEFINE VARIABLE-----------------------------------------------------------------%%
    % Flow Controls
    global AOA_EST_MODE SIMULATION
    
    % Debug Controls
    global NUMBER_OF_PACKETS_TO_CONSIDER
    
    % Physical parameter
    global n_antenna n_subcarrier
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %% Packet reshaping, sanitize and sampling
    if ~SIMULATION
        csi_trace = cell(NUMBER_OF_PACKETS_TO_CONSIDER,1);
        for i = 1:NUMBER_OF_PACKETS_TO_CONSIDER
            csi_trace{i}.csi = zeros(n_antenna,n_subcarrier);
        end
        for i = 1:n_antenna
            tmp_trace = read_log_file([data_file '_' num2str(i)]);
            for j = 1:NUMBER_OF_PACKETS_TO_CONSIDER
                csi_trace{j}.csi(i,:) = tmp_trace{j}.csi(1,1,:);
            end
        end
    else
        load(data_file);
    end

    % Set the number of packets to consider, by default consider all
    num_packets = length(csi_trace);
    if NUMBER_OF_PACKETS_TO_CONSIDER ~= -1
        num_packets = NUMBER_OF_PACKETS_TO_CONSIDER;
    end

    sampled_csi_trace = csi_sampling(csi_trace, num_packets, ...
            1, length(csi_trace));
    sanitized_csi = cell(NUMBER_OF_PACKETS_TO_CONSIDER,1);
    
    csi1 = sampled_csi_trace{1}.csi;
    [sanitized_csi{1}.csi, est_slope] = spotfi_algorithm_1(csi1);
    
    for i = 2:NUMBER_OF_PACKETS_TO_CONSIDER
        tmp_csi = sampled_csi_trace{i}.csi;
        sanitized_csi{i}.csi = spotfi_algorithm_1(tmp_csi, est_slope);
    end   
    
    %% mode choosing
    if strcmp(AOA_EST_MODE, 'SPOTFI')
        get_aoa = @(csi_packet)spotfi(csi_packet.csi);
    elseif strcmp(AOA_EST_MODE, 'MUSIC')
        get_aoa = @(csi_packet, ind_subcarrier)musicAOA(csi_packet.csi(:,ind_subcarrier));
    else
        error('No available algorithm specify')
    end
    
    %% enter your main function here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_top_aoa = get_aoa(sanitized_csi{2},3);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end