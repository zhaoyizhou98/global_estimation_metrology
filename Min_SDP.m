function varargout = Min_SDP(strategy,dimensions,Phi,H)
% MIN_SDP 
% strategy: 1 for i, 2 for ii, 3 for iii, 4 for iv
% dimensions = [I1,O1,I2,O2]
% 
% result is of form [obj,tildeY,lambda,h] for i, ii and iv
% result is of form [obj,Y1,Y2,lambda,h] for iii
% 
% We only consider using two signal channels for iii and iv 
% and using two or three channels for i and ii. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
shape = size(Phi);
dim = shape(1);
components = shape(2);
N = length(dimensions)/2;

switch strategy
    % 
    % Strategies of type i
    % 
    case 1
        cvx_begin sdp quiet
            % variables
            variable tildeY(dim,dim) hermitian;
            variable lambda;
            variable h(components,components) hermitian;
            expression bigmat(dim+components,dim+components);
            dimout = prod(dimensions(2*(1:N)));
            bigmat(1:dim,1:dim) = lambda*eye(dim)/dimout+tildeY;
            bigmat(1:dim,dim+1:end) = 2*(conj(H)*conj(Phi)-1j*conj(Phi)*h);
            bigmat(dim+1:end,1:dim) = bigmat(1:dim,dim+1:end)';
            bigmat(dim+1:end,dim+1:end) = eye(components);
            % restrictions
            bigmat >= 0;
            LambdaY = myoperation(tildeY,2*(1:N),dimensions);
            LambdaY == 0;
            % objective function
            minimize(lambda);
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = tildeY;
        varargout{3} = lambda;
        varargout{4} = h;
    % 
    % Strategies of type ii
    % 
    case 2
        cvx_begin sdp quiet
            % variables
            variable tildeY(dim,dim) hermitian;
            variable lambda;
            variable h(components,components) hermitian;
            expression bigmat(dim+components,dim+components);
            dimout = prod(dimensions(2*(1:N)));
            bigmat(1:dim,1:dim) = lambda*eye(dim)/dimout+tildeY;
            bigmat(1:dim,dim+1:end) = 2*(conj(H)*conj(Phi)-1j*conj(Phi)*h);
            bigmat(dim+1:end,1:dim) = bigmat(1:dim,dim+1:end)';
            bigmat(dim+1:end,dim+1:end) = eye(components);
            % restrictions
            bigmat >= 0;
            if N == 2
                o2Y = myoperation(tildeY,4,dimensions);
                i2o2Y = myoperation(tildeY,[3,4],dimensions);
                o1i2o2Y = myoperation(tildeY,[2,3,4],dimensions);
                o2Y - i2o2Y +o1i2o2Y == 0;
            elseif N == 3
                o3Y = myoperation(tildeY,6,dimensions);
                i3o3Y = myoperation(tildeY,[5,6],dimensions);
                o2i3o3Y = myoperation(tildeY,[4,5,6],dimensions);
                i2o2i3o3Y = myoperation(tildeY,[3,4,5,6],dimensions);
                o1i2o2i3o3Y = myoperation(tildeY,[2,3,4,5,6],dimensions);
                o3Y - i3o3Y + o2i3o3Y - i2o2i3o3Y + o1i2o2i3o3Y == 0;
            else
                error("Only support N = 2 or N = 3.");
            end
            % objective function
            minimize(lambda);
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = tildeY;
        varargout{3} = lambda;
        varargout{4} = h;
    % 
    % Strategies of type iii
    % 
    case 3
        cvx_begin sdp quiet
        % variables
        variable Y1(dim,dim) hermitian;
        variable Y2(dim,dim) hermitian;
        variable lambda;
        variable h(components,components) hermitian;
        expression bigmat1(dim+components,dim+components);
        expression bigmat2(dim+components,dim+components);
        dimout = prod(dimensions(2*(1:N)));
        bigmat1(1:dim,1:dim) = lambda*eye(dim)/dimout+Y1;
        bigmat1(1:dim,dim+1:end) = 2*(conj(H)*conj(Phi)-1j*conj(Phi)*h);
        bigmat1(dim+1:end,1:dim) = bigmat1(1:dim,dim+1:end)';
        bigmat1(dim+1:end,dim+1:end) = eye(components);

        bigmat2(1:dim,1:dim) = lambda*eye(dim)/dimout+Y2;
        bigmat2(1:dim,dim+1:end) = 2*(conj(H)*conj(Phi)-1j*conj(Phi)*h);
        bigmat2(dim+1:end,1:dim) = bigmat2(1:dim,dim+1:end)';
        bigmat2(dim+1:end,dim+1:end) = eye(components);
        % restrictions
        bigmat1 >= 0; bigmat2 >= 0;
        o2Y1 = myoperation(Y1,4,dimensions);
        i2o2Y1 = myoperation(Y1,[3,4],dimensions);
        o1i2o2Y1 = myoperation(Y1,[2,3,4],dimensions);
        o2Y1 - i2o2Y1 +o1i2o2Y1 == 0;

        o1Y2 = myoperation(Y2,2,dimensions);
        i1o1Y2 = myoperation(Y2,[1,2],dimensions);
        i1o1o2Y2 = myoperation(Y2,[1,2,4],dimensions);
        o1Y2 - i1o1Y2+ i1o1o2Y2 == 0;
        % objective function
        minimize(lambda);
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = Y1;
        varargout{3} = Y2;
        varargout{4} = lambda;
        varargout{5} = h;
    % 
    % Strategies of type iv
    % 
    case 4
        cvx_begin sdp quiet
            % variables
            variable tildeY(dim,dim) hermitian;
            variable lambda;
            variable h(components,components) hermitian;
            expression bigmat(dim+components,dim+components);
            dimout = prod(dimensions(2*(1:N)));
            bigmat(1:dim,1:dim) = lambda*eye(dim)/dimout+tildeY;
            bigmat(1:dim,dim+1:end) = 2*(conj(H)*conj(Phi)-1j*conj(Phi)*h);
            bigmat(dim+1:end,1:dim) = bigmat(1:dim,dim+1:end)';
            bigmat(dim+1:end,dim+1:end) = eye(components);
            % restrictions
            bigmat >= 0;
            i1o1o2Y = myoperation(tildeY,[1,2,4],dimensions);
            i1o1Y = myoperation(tildeY,[1,2],dimensions);
            o1i2o2Y = myoperation(tildeY,[2,3,4],dimensions);
            i2o2Y = myoperation(tildeY,[3,4],dimensions);
            o1Y = myoperation(tildeY,2,dimensions);
            o2Y = myoperation(tildeY,4,dimensions);
            o1o2Y = myoperation(tildeY,[2,4],dimensions);
            i1o1o2Y- i1o1Y +o1i2o2Y -i2o2Y + o1Y+o2Y-o1o2Y == 0;
            % objective function
            minimize(lambda);
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = tildeY;
        varargout{3} = lambda;
        varargout{4} = h;
    otherwise
        error("Input strategy is illegal. 1 for i, 2 for ii, 3 for iii, 4 for iv.");
end