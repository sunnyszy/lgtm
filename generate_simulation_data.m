%generate the simulation data
function generate_simulation_data(file_name, snr, n_packet)

if nargin < 3
    n_packet = 10;
end

global d theta_the l_the channel_frequency delta_f n_subcarrier c;

base_frequency = channel_frequency - delta_f * (n_subcarrier/2 - 0.5);
x_offsets = [d 0 -d];
tracking_x = l_the*cos(theta_the);
tracking_y = l_the*sin(theta_the);


%add multipath
wall1_x = -8;
wall2_x = 8;

%add movement
v = 0.01; %0.01m/s


csi_trace = cell(n_packet,1);

%generate each pkt
for p = 1:n_packet
    tracking_x = tracking_x + 2*v*rand()-0.01;
    tracking_y = tracking_y + 2*v*rand()-0.01;
    
    mirror_x1 = 2*wall1_x - tracking_x;
    mirror_x2 = 2*wall2_x - tracking_x;
    mirror_y1 = tracking_y;
    mirror_y2 = tracking_y;
    csi = zeros(3,n_subcarrier);
    for i = 1:3
        x_offset = x_offsets(i);
        for j = 1:n_subcarrier
            tmp_frequency = base_frequency + delta_f * (j-1);
            distance = norm([tracking_x-x_offset tracking_y]);
            distance_mir1 = norm([mirror_x1-x_offset mirror_y1]);
            distance_mir2 = norm([mirror_x2-x_offset mirror_y2]);
            csi(i, j) = (1/distance)*exp(-1i*2*pi*tmp_frequency/c*distance)...
                + (1/distance_mir1)*exp(-1i*2*pi*tmp_frequency/c*distance_mir1) ...
                + (1/distance_mir2)*exp(-1i*2*pi*tmp_frequency/c*distance_mir2);
        end
    end

    %add noise to generate n packet
    for i = 1:n_packet
        csi_trace{i}.csi = zeros(1,3,n_subcarrier);
        csi_trace{i}.csi(1,:,:)= awgn(csi,snr, 'measured');
    end
end
save(file_name, 'csi_trace');    

end