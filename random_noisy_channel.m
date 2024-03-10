% Randomly generate noise channels to check the strict hierarchy between
% different types of strategies

clear all; clc;
% number of channels to be generated
rep_len = 1000;
obj1 = zeros(rep_len,1);
obj2 = zeros(rep_len,1);
obj3 = zeros(rep_len,1);
obj4 = zeros(rep_len,1);

syms theta real;
N = 2;
int_low = -pi;
int_high = pi;
prior = 1/(2*pi);
% tolerance
tol = 1e-04;

% parallel for
pool = parpool(8);
parfor i = 1:rep_len
    [Phi,H] = get(prior,N,int_low,int_high);
    obj1(i) = Min_SDP(1,2*ones(1,2*N),Phi,H);
    obj2(i) = Min_SDP(2,2*ones(1,2*N),Phi,H);
    obj3(i) = Min_SDP(3,2*ones(1,2*N),Phi,H);
    obj4(i) = Min_SDP(4,2*ones(1,2*N),Phi,H);
end
delete(pool);
save(sprintf('./Data/random_noise_N%d.mat',N),'obj1','obj2','obj3','obj4');

function [Phi, H] = get(prior,N,int_low,int_high)
    % in fact, N can only be 2, since we mainly care about N = 2 case
    if N ~= 2
        error("N must be 2")
    end
    syms theta real;
    U = [exp(-1j*theta/2) 0;0 exp(1j*theta/2)];
    Kraus = myrandomChannel(2);
    K_num = 2;
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