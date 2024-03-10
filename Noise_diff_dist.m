% Noise channel with different prior distributions
clear all; clc;
syms theta real;
N = 2;
int_low = -pi; int_high = pi;
eta = 0.5;
gamma = 0.7;
rep_len = 100;
obj1 = zeros(rep_len,1);
obj2 = zeros(rep_len,1);
obj3 = zeros(rep_len,1);
obj4 = zeros(rep_len,1);


% Gaussian distribution with different standard deviation [0.1,2]
% parallel for
pool = parpool(8);
step = (2-0.1)/rep_len;
parfor i = 1:rep_len
    delta = 0.1+i*step;
    prior = exp(-theta^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    [Phi,H] = get(prior,N,int_low,int_high,gamma,eta);
    obj1(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
    obj2(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    obj3(i) = Min_SDP(3,2*ones(1,2*N),Phi,H);
    obj4(i) = Min_SDP(4,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/noise_gaussian_diffdelta_N%d.mat',N),'obj1','obj2','obj3','obj4');

% Gaussian distribution with different mean value [-3,3]
pool = parpool(8);
delta = 2;
step = 6/rep_len;
parfor i = 1:rep_len
    mu = -3+i*step;
    prior = exp(-(theta-mu)^2/(2*delta^2));
    prior = prior/double(int(prior,theta,int_low,int_high));
    [Phi,H] = get(prior,N,int_low,int_high,gamma,eta);
    obj1(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
    obj2(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    obj3(i) = Min_SDP(3,2*ones(1,2*N),Phi,H);
    obj4(i) = Min_SDP(4,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/noise_gaussian_diffmean_N%d.mat',N),'obj1','obj2','obj3','obj4');

% Gaussian mixture distribution with different weight
pool = parpool(8);
delta1 = 1; mu1 = -pi/2;
delta2 = 2; mu2 = pi/2;
step = 1/rep_len;
parfor i = 1:rep_len
    w = i*step;
    prior = w*exp(-(theta-mu1)^2/(2*delta1^2))/sqrt(2*pi*delta1^2)+(1-w)*exp(-(theta-mu2)^2/(2*delta2^2))/sqrt(2*pi*delta2^2);
    prior = prior/double(int(prior,theta,int_low,int_high));
    [Phi,H] = get(prior,N,int_low,int_high,gamma,eta);
    obj1(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
    obj2(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    obj3(i) = Min_SDP(3,2*ones(1,2*N),Phi,H);
    obj4(i) = Min_SDP(4,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/noise_GM_diffw_N%d.mat',N),'obj1','obj2','obj3','obj4');

% Beta distribution with different a
pool = parpool(8);
b = 2;
step = (2-0.2)/rep_len;
parfor i = 1:rep_len
    a = 0.2+i*step;
    prior = ((theta+pi)/(2*pi))^(a-1)*(1-(theta+pi)/(2*pi))^(b-1);
    prior = prior/integral(matlabFunction(prior),int_low,int_high);
    [Phi,H] = get2(prior,N,int_low,int_high,gamma,eta);
    obj1(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
    obj2(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    obj3(i) = Min_SDP(3,2*ones(1,2*N),Phi,H);
    obj4(i) = Min_SDP(4,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/noise_beta_diffa_N%d.mat',N),'obj1','obj2','obj3','obj4');


function [Phi, H] = get(prior,N,int_low,int_high,gamma,eta)
    % in fact, N can only be 2, since we mainly care about N = 2 case
    if N ~= 2
        error("N must be 2")
    end
    syms theta real;
    U = [exp(-1j*theta/2) 0;0 exp(1j*theta/2)];
    Kraus(:,:,1) = [1 0;0 sqrt(1-gamma)]*sqrt(eta)*eye(2);
    Kraus(:,:,2) = [0 sqrt(gamma);0 0]*sqrt(eta)*eye(2);
    Kraus(:,:,3) = [1 0;0 sqrt(1-gamma)]*sqrt(1-eta)*[0 1;1 0];
    Kraus(:,:,4) = [0 sqrt(gamma);0 0]*sqrt(1-eta)*[0 1;1 0];
    K_num = 4;
    K = sym(zeros(2,2,K_num));
    for i = 1:K_num
        K(:,:,i) = Kraus(:,:,i)*U;
    end
    % create Phi and H
    ctheta = zeros(4^N,4^N);
    for i = 1:K_num
        for j = 1:K_num
            ctheta = ctheta + kron(reshape(K(:,:,i),[4,1]), reshape(K(:,:,j),[4,1]))* ... 
                kron(reshape(K(:,:,i),[4,1]), reshape(K(:,:,j),[4,1]))';
        end
    end
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

% a more efficient numerical integral
function [Phi, H] = get2(prior,N,int_low,int_high,gamma,eta)
    % in fact, N can only be 2, since we mainly care about N = 2 case
    if N ~= 2
        error("N must be 2")
    end
    syms theta real;
    U = [exp(-1j*theta/2) 0;0 exp(1j*theta/2)];
    Kraus(:,:,1) = [1 0;0 sqrt(1-gamma)]*sqrt(eta)*eye(2);
    Kraus(:,:,2) = [0 sqrt(gamma);0 0]*sqrt(eta)*eye(2);
    Kraus(:,:,3) = [1 0;0 sqrt(1-gamma)]*sqrt(1-eta)*[0 1;1 0];
    Kraus(:,:,4) = [0 sqrt(gamma);0 0]*sqrt(1-eta)*[0 1;1 0];
    K_num = 4;
    K = sym(zeros(2,2,K_num));
    for i = 1:K_num
        K(:,:,i) = Kraus(:,:,i)*U;
    end
    % create Phi and H
    ctheta = zeros(4^N,4^N);
    for i = 1:K_num
        for j = 1:K_num
            ctheta = ctheta + kron(reshape(K(:,:,i),[4,1]), reshape(K(:,:,j),[4,1]))* ... 
                kron(reshape(K(:,:,i),[4,1]), reshape(K(:,:,j),[4,1]))';
        end
    end
    cbar = integral(matlabFunction(ctheta*prior),int_low,int_high,'ArrayValued',true);
    thetacbar = integral(matlabFunction(ctheta*prior*theta),int_low,int_high,'ArrayValued',true);
    H = mysylvester(-cbar,-cbar,thetacbar);
    [V D] = eig(cbar);
    zero_idx = find(abs(diag(D)) <= 0.00001);
    Phi = V*sqrt(D);
    Phi(:,zero_idx) = [];
end