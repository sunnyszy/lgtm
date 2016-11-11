clc,clear;


%slow walk
csi_profile = read_log_file('./test-data/d2.dat');
n_entry = length(csi_profile);
tmp_rssi1 = zeros(1, n_entry);
for j = 1:n_entry
    tmp_rssi1(j) = csi_profile{j}.rssi3;
end
plot(tmp_rssi1, 'b');
xlabel('sequence');
ylabel('RSSI/DB?');
% xlim([0,1500]);
