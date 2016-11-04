%generate the simulation data
function generate_simulation_data(file_name, snr, n_packet)

if nargin == 2
    n_packet = 10;
end


global d theta_the l_the channel_frequency delta_f n_subcarrier c;

base_frequency = channel_frequency - delta_f * (n_subcarrier/2 - 0.5);
csi = zeros(1,3, n_subcarrier);
x_offsets = [d 0 -d];

%generate one packet
for i = 1:3
    x_offset = x_offsets(i);
    for j = 1:n_subcarrier
        tmp_frequency = base_frequency + delta_f * (j-1);
        distance = norm([l_the*cos(theta_the)-x_offset l_the*sin(theta_the)]);
        csi(1, i, j) = (1/distance)*exp(-1i*2*pi*tmp_frequency/c*distance);
    end
end

%add noise to generate n packet
squeezed_csi = squeeze(csi);
csi_trace = cell(n_packet,1);
for i = 1:n_packet
    csi_trace{i}.csi = zeros(1,3,n_subcarrier);
    csi_trace{i}.csi(1,:,:)= awgn(squeezed_csi,snr, 'measured');
end

save(file_name, 'csi_trace');    

end