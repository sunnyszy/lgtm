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

%% Time of Flight (ToF) Sanitization Algorithm, find a linear fit for the unwrapped CSI phase
% csi_matrix -- the CSI matrix whose phase is to be adjusted
% delta_f    -- the difference in frequency between subcarriers
% Return:
% csi_matrix -- the same CSI matrix with modified phase
function [csi_matrix, est_slope] = spotfi_algorithm_1(csi_matrix, est_slope)
    %% Time of Flight (ToF) Sanitization Algorithm
    % Unwrap phase from CSI matrix
    global n_subcarrier n_antenna
    R = abs(csi_matrix);
    phase_matrix = unwrap(angle(csi_matrix), pi, 2);

    % Parse input args, if < 2: the first packet, fit STO
    if nargin < 2
%         fit_X = zeros(n_antenna*n_subcarrier, 1);
%         fit_Y = zeros(n_antenna*n_subcarrier, 1);
        tmp_X = repmat(1:n_subcarrier,[n_antenna,1]); % correspond X
        fit_X = reshape(tmp_X, 1, n_subcarrier*n_antenna);
        fit_Y = reshape(phase_matrix,1, n_antenna*n_subcarrier);

        % Linear fit is common across all antennas
        para = polyfit(fit_X, fit_Y, 1);
        est_slope = para(1);
    end
    
    % correspond subtraction slope
    tmp_X = repmat(1:n_subcarrier,[n_antenna,1]);
    tmp_X = tmp_X - 1;
    phase_matrix = phase_matrix - tmp_X*est_slope;
    
    
    % Reconstruct the CSI matrix with the adjusted phase
    % TODO: for packet 2-m, only change R?
    csi_matrix = R .* exp(1i * phase_matrix);
end