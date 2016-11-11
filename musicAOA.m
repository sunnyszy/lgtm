%estimate theta by csi information, the labels are one a line
%h: csi information, M*1 vector, ith item is the csi of ith antenna
%l: the distance of two near label; lambda: lambda of the wave
%M: number of labels; D: number of eigenvector to throw away
%P: a 1*180 vector -- ith item represent the possibility of the angle in i
%degree (unnormalized)
function [top_aoa] = musicAOA(h)

global c channel_frequency d n_antenna

lambda = c/channel_frequency;

R = h*h'; 
% Find the eigenvalues and eigenvectors of the covariance matrix
[eigenvectors, eigenvalue_matrix] = eig(R);

if ~issorted(diag(eigenvalue_matrix))
    [eigenvalue_matrix,I] = sort(diag(eigenvalue_matrix));
    eigenvectors = eigenvectors(:, I);
end

% Find max eigenvalue for normalization
max_eigenvalue = eigenvalue_matrix(n_antenna,n_antenna);

for ii = 1:size(n_antenna, 1)
    eigenvalue_matrix(ii, ii) = eigenvalue_matrix(ii, ii) / max_eigenvalue;
end

% Find the largest decrease ratio that occurs
ratio = zeros(1, n_antenna -2);
for i = 1:n_antenna-2
    if eigenvalue_matrix(i,i) == 0
        ratio(i) = 0;
    else
        ratio(i) = eigenvalue_matrix(i+1,i+1)/eigenvalue_matrix(i,i);
    end
end

[~, ind] = max(ratio);
%%


E = eigenvectors(:, 1:n_antenna - ind);
P = zeros(1, 180);
E_product = E*E';


for theta=1:180
    steering_vector = exp(-1i*((n_antenna:-1:1)-1)*2*pi*d/lambda*cos(deg2rad(theta)));
    P(theta)=abs(1./(steering_vector*E_product*(steering_vector)'));
end

[~,top_aoa] = max(P);

end

