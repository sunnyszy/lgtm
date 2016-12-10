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

function wgtt_runner()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    globals_init()
    
    %flow control
    global SIMULATION SIMULAIION_ALWAYS_GENERATE_DATA AOA_EST_MODE OUTDOOR_FLAG
    SIMULATION = false;
    SIMULAIION_ALWAYS_GENERATE_DATA = true;
    AOA_EST_MODE = 'SPOTFI';%'MUSIC' 'SPOTFI';
    OUTDOOR_FLAG = false;
    
    %output control
    global OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH 
    OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH  = true;
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
        name_base = '90';
    end

    data_file = ['test-data/' name_base];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% mode choosing
    % read data algorithm
    if ~SIMULATION && OUTDOOR_FLAG %outdoor experiment
        load_csi_trace = @(data_file) outdoor_load_csi_trace(data_file);
        main = @(csi_trace,get_aoa) outdoor_main(csi_trace,get_aoa);
    elseif  ~SIMULATION && ~OUTDOOR_FLAG % indoor experiment
        load_csi_trace = @(data_file) indoor_load_csi_trace(data_file);
        main = @(csi_trace,get_aoa) indoor_main(csi_trace,get_aoa);
    elseif SIMULATION
        load_csi_trace = @(data_file) load(data_file);
        main = @(csi_trace,get_aoa) simulation_main(csi_trace,get_aoa);
    end
    % aoa algorithm
    if strcmp(AOA_EST_MODE, 'SPOTFI')
        get_aoa = @(csi_packet)spotfi(csi_packet.csi);
    elseif strcmp(AOA_EST_MODE, 'MUSIC')
        get_aoa = @(csi_packet, ind_subcarrier)musicAOA(csi_packet.csi(:,ind_subcarrier));
    else
        error('No available algorithm specify')
    end
    
    %% start running
    % load file
    csi_trace = load_csi_trace(data_file);
    
    %% Packet reshaping, sanitize and sampling
    % Set the number of packets to consider, by default consider all
    num_packets = length(csi_trace);
    sanitized_csi = cell(num_packets,1);
    
    csi1 = csi_trace{1}.csi;
    [sanitized_csi{1}.csi, est_slope] = spotfi_algorithm_1(csi1);
    
    for i = 2:num_packets
        tmp_csi = csi_trace{i}.csi;
        sanitized_csi{i}.csi = spotfi_algorithm_1(tmp_csi, est_slope);
    end   
    
    
    
    % main function here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    main(sanitized_csi,get_aoa);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(['test-output/' name_base]);
end

%% csi trace read for outdoor
function csi_trace = outdoor_load_csi_trace(data_file_name)
% Physical parameter
global n_antenna


%% enter your read data function here
tmp_trace = read_log_file([data_file_name '.dat']);
n_pkts = length(tmp_trace);
% csi_trace = cell(n_pkts,1);
counter = 1;
% for i = 1:n_pkts
%     csi_trace{i}.csi = zeros(n_antenna,n_subcarrier);
% end
for i = 1:n_pkts
    if tmp_trace{i}.csi_len 
%         tmp_trace{i}.csi_len
        csi_noempty{counter,1}.csi = squeeze(tmp_trace{i}.csi(3,1,:));
        counter = counter + 1;
    end
end

start_index = 16950;
end_index = 17050;
range = end_index - start_index + 1;
n_aggre = range - n_antenna + 1;
csi_trace = cell(n_aggre, 1);
for i = 1:n_aggre
    csi_trace{i}.csi = [];
    for j = 1:n_antenna
        csi_trace{i}.csi = [csi_trace{i}.csi; transpose(csi_noempty{i+start_index+j-2}.csi)];
    end
end

end


%% csi_trace read for indoor test
function csi_trace = indoor_load_csi_trace(data_file_name)

% Physical parameter
global n_antenna
max_n_pkt = -1;
%% enter your read data function here
for i = 1:n_antenna
    tmp_trace = read_log_file([data_file_name '_' num2str(i)]);
    if length(tmp_trace) > max_n_pkt
        max_n_pkt = length(tmp_trace);
    end
end

%load all pkts (possibably with no csi)
noempty_trace = cell(n_antenna, max_n_pkt);
min_n_pkt = 1e8;
for i = 1:n_antenna
    tmp_trace = read_log_file([data_file_name '_' num2str(i)]);
    counter = 1;
    for j = 1:length(tmp_trace)
        if tmp_trace{j}.csi_len
            noempty_trace{i,counter} = tmp_trace{j}.csi(1,1,:);
            counter = counter + 1;
        end
    end
    if (counter-1) < min_n_pkt
        min_n_pkt = counter - 1;
    end
end

csi_trace = cell(min_n_pkt, 1);
for i = 1:min_n_pkt
    for j = 1:n_antenna
        csi_trace{i}.csi(j,:) = noempty_trace{j,i};
    end
end

end

%% indoor main
function indoor_main(csi_trace, get_aoa)
    n_pkt = 10;
    aoas = zeros(n_pkt, 1);
    for i = 1:n_pkt
        aoas(i) = get_aoa(csi_trace{i});
        disp(i);
    end
    plot(aoas);
    legend(['mean aoa: ' num2str(mean(aoas))]);
end

%% outdoor main
function outdoor_main(csi_trace, get_aoa)
    global d;
    d = 0.01;
    n_pkt = 10;
    aoas = zeros(n_pkt, 1);
    for i = 1:n_pkt
        aoas(i) = get_aoa(csi_trace{i});
        disp(i);
    end
    plot(aoas);
end
