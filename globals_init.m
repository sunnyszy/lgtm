function globals_init
    %% DEBUG AND OUTPUT VARIABLES-----------------------------------------------------------------%%
    % Flow Controls
    global AOA_EST_MODE SIMULATION SIMULAIION_ALWAYS_GENERATE_DATA ...
        OUTDOOR_FLAG
    AOA_EST_MODE = 'SPOTFI'; % 'MUSIC'  'SPOTFI'
    SIMULATION = true;
    SIMULAIION_ALWAYS_GENERATE_DATA = false;
    OUTDOOR_FLAG = false;
%     NUMBER_OF_PACKETS_TO_CONSIDER = 10; % Set to -1 to ignore this variable's value    
    
    % Output controls
    global OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH
    OUTPUT_AOA_TOF_MUSIC_PEAK_GRAPH = false;
    
    % constant parameter
    global d theta_the l_the channel_frequency delta_f n_subcarrier c ...
    n_antenna subset_x subset_y;
    theta_the = pi/6; % in rad
    l_the = 5; % in meter
    channel_frequency = 2462e6;
    delta_f = 312.5e3;
    n_subcarrier = 56; % Atheros 56
    c = 3e8; % speed of light
    n_antenna = 3;
    d = c/channel_frequency/2; % distance between two antenna;
    [subset_x, subset_y] = best_rank_subset(); %antenna subset used in smooth matrix
    
    function [sub_x, sub_y] = best_rank_subset()
        score = zeros(n_subcarrier, n_antenna);
        for x = 1:n_subcarrier-1
            for y = 1:n_antenna-1
                if x*y <= (n_antenna+1-y)*(n_subcarrier+1-x) % for the best of full rank
                    score(x,y) = x*y;
                else
                    score(x,y) = 0;
                end
            end
        end
        ind = find(score==max(score(:)));
        [sub_x, sub_y] = ind2sub([n_subcarrier, n_antenna],ind);
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end