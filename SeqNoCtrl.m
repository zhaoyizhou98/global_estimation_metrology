% Unexpected result N = 2
% Sequential strategies without control
clear all; clc;
int_low = -pi; int_high = pi;
prior = 1/(2*pi);
[Phi, H] = get(prior, int_low, int_high);
% without control is equivalent to parallel strategies
obj = Min_SDP(1,2*ones(1,2),Phi,H);
obj

% get function for without control N = 2
function [Phi, H] = get(prior,int_low,int_high)
    syms theta real;
    U = [exp(-1j*theta) 0;0 exp(1j*theta)];
    vecU = reshape(U,[4,1]);
    ctheta = vecU*vecU';
    % create Phi and H
    cbar = int(ctheta*prior,theta,int_low,int_high);
    thetacbar = int(ctheta*prior*theta,theta,int_low,int_high);
    cbar = double(cbar);
    thetacbar = double(thetacbar);
    H = mysylvester(-cbar,-cbar,thetacbar);
    [V D] = eig(cbar);
    zero_idx = find(abs(diag(D)) <= 0.00001);
    Phi = V*sqrt(D);
    Phi(:,zero_idx) = [];
end