% Impact of different prior distributions
clear all; clc;
syms theta real;
int_low = -pi; int_high = pi;
rep_len = 100;
objglobal = zeros(rep_len,1);
objlocal = zeros(rep_len,1);
step = (2-0.1)/rep_len;

% N = 2
N = 2;
rho = sym(zeros(4,4));
rho(1,1) = 1/2; rho(1,4) = exp(-2j*theta)/2; rho(4,1) = exp(2j*theta)/2; rho(4,4) = 1/2;
pool = parpool(8);
parfor i = 1:rep_len
    delta = step*i+0.1;
    prior = exp(-theta^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    rhobar = int(rho*prior,theta,int_low,int_high);
    thetarhobar = int(rho*theta*prior,theta,int_low,int_high);
    rhobar = double(rhobar); thetarhobar = double(thetarhobar);
    S = mysylvester(rhobar,rhobar,2*thetarhobar);
    objlocal(i) = trace(rhobar*S*S);
    [Phi, H] = get(prior,N,int_low,int_high);
    objglobal(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/GHZ_diffdelta_N%d.mat',N),'objglobal','objlocal');

% N = 3
N = 3;
rho = sym(zeros(8,8));
rho(1,1) = 1/2; rho(1,8) = exp(-3j*theta)/2; rho(8,1) = exp(3j*theta)/2; rho(8,8) = 1/2;
pool = parpool(8);
parfor i = 1:rep_len
    delta = step*i+0.1;
    prior = exp(-theta^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    rhobar = int(rho*prior,theta,int_low,int_high);
    thetarhobar = int(rho*theta*prior,theta,int_low,int_high);
    rhobar = double(rhobar); thetarhobar = double(thetarhobar);
    S = mysylvester(rhobar,rhobar,2*thetarhobar);
    objlocal(i) = trace(rhobar*S*S);
    [Phi, H] = get(prior,N,int_low,int_high);
    objglobal(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/GHZ_diffdelta_N%d.mat',N),'objglobal','objlocal');

% With or without control N = 2
N = 2;
pool = parpool(8);
parfor i = 1:rep_len
    delta = step*i+0.1;
    prior = exp(-theta^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    [Phi, H] = get(prior, N, int_low, int_high);
    objglobal(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    [Phi, H] = get2(prior, int_low, int_high);
    objlocal(i) = Min_SDP(1,2*ones(1,2),Phi,H);
end
delete(pool);
save(sprintf('./Data/SeqCtrl_diffprior_N%d.mat',N),'objglobal','objlocal');

% With or without control N = 3
N = 3;
pool = parpool(8);
parfor i = 1:rep_len
    delta = step*i+0.1;
    prior = exp(-theta^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    [Phi, H] = get(prior, N, int_low, int_high);
    objglobal(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    [Phi, H] = get3(prior, int_low, int_high);
    objlocal(i) = Min_SDP(1,2*ones(1,2),Phi,H);
end
delete(pool);
save(sprintf('./Data/SeqCtrl_diffprior_N%d.mat',N),'objglobal','objlocal');


function [Phi, H] = get(prior,N,int_low,int_high)
    syms theta real;
    U = [exp(-1j*theta/2) 0;0 exp(1j*theta/2)];
    vecU = reshape(U,[4,1]);
    if N == 2
        ctheta = kron(vecU,vecU);
    elseif N == 3
        ctheta = kron(vecU,kron(vecU,vecU));
    else
        error("N must be 2 or 3");
    end
    ctheta = ctheta*ctheta';
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
% get function for without control N = 2
function [Phi, H] = get2(prior,int_low,int_high)
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
% get function for without control N = 3
function [Phi, H] = get3(prior,int_low,int_high)
    syms theta real;
    U = [exp(-3j*theta/2) 0;0 exp(3j*theta/2)];
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