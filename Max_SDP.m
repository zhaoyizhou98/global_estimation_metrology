function varargout = Max_SDP(strategy,dimensions,Phi,H)
% MAX_SDP 
% strategy: 1 for i, 2 for ii, 3 for iii, 4 for iv
% dimensions = [I1,O1,I2,O2]
% 
% result is of form [obj,X,B,C] for i, ii and iv
% result is of form [obj,X1,X2,B,C] for iii
% 
% We only consider using two signal channels for iii and iv
% and using two or three channels for i and ii.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
            variable X(dim,dim) complex semidefinite;
            variable C(components,components) complex semidefinite;
            variable B(components,dim) complex;
            expression bigmat(components+dim,components+dim);
            bigmat(1:dim,1:dim) = X;
            bigmat(dim+1:end,1:dim) = B;
            bigmat(1:dim,dim+1:end) = B';
            bigmat(dim+1:end,dim+1:end) = C;
            % restrictions
            bigmat >= 0;
            B*conj(Phi) == Phi.'*B';
            oX = myoperation(X,2*(1:N),dimensions);
            X == oX;
            trace(X) == prod(dimensions(2*(1:N)));
            % objective function
            maximize(-trace(C)-4*real(trace(conj(Phi)*B*conj(H))));
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = X;
        varargout{3} = B;
        varargout{4} = C;
    % 
    % Strategies of type ii
    % 
    case 2
        cvx_begin sdp quiet
            % variables
            variable X(dim,dim) complex semidefinite;
            variable C(components,components) complex semidefinite;
            variable B(components,dim) complex;
            expression bigmat(components+dim,components+dim);
            bigmat(1:dim,1:dim) = X;
            bigmat(dim+1:end,1:dim) = B;
            bigmat(1:dim,dim+1:end) = B';
            bigmat(dim+1:end,dim+1:end) = C;
            % restrictions
            bigmat >= 0;
            B*conj(Phi) == Phi.'*B';
            trace(X) == prod(dimensions(2*(1:N)));
            if N == 2
                o2X1 = myoperation(X,4,dimensions);
                i2o2X1 = myoperation(X,[3,4],dimensions);
                o1i2o2X1 = myoperation(X,[2,3,4],dimensions);
                X == o2X1 - i2o2X1 +o1i2o2X1;
            elseif N == 3
                o3X = myoperation(X,6,dimensions);
                i3o3X = myoperation(X,[5,6],dimensions);
                o2i3o3X = myoperation(X,[4,5,6],dimensions);
                i2o2i3o3X = myoperation(X,[3,4,5,6],dimensions);
                o1i2o2i3o3X = myoperation(X,[2,3,4,5,6],dimensions);
                X == o3X - i3o3X + o2i3o3X - i2o2i3o3X + o1i2o2i3o3X;
            else
                error("Only support N = 2 or N = 3.");
            end
            % objective function
            maximize(-trace(C)-4*real(trace(conj(Phi)*B*conj(H))));
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = X;
        varargout{3} = B;
        varargout{4} = C;
    % 
    % Strategies of type iii
    % 
    case 3
        cvx_begin sdp quiet
            % variables
            variable X1(dim,dim) complex semidefinite;
            variable X2(dim,dim) complex semidefinite;
            variable C(components,components) complex semidefinite;
            variable B(components,dim) complex;
            expression bigmat(components+dim,components+dim);
            bigmat(1:dim,1:dim) = X1+X2;
            bigmat(dim+1:end,1:dim) = B;
            bigmat(1:dim,dim+1:end) = B';
            bigmat(dim+1:end,dim+1:end) = C;
            % restrictions
            bigmat >= 0;
            B*conj(Phi) == Phi.'*B';
            trace(X1) + trace(X2) == dimensions(2)*dimensions(4);
            o2X1 = myoperation(X1,4,dimensions);
            i2o2X1 = myoperation(X1,[3,4],dimensions);
            o1i2o2X1 = myoperation(X1,[2,3,4],dimensions);
            X1 == o2X1 - i2o2X1 +o1i2o2X1;

            o1X2 = myoperation(X2,2,dimensions);
            i1o1X2 = myoperation(X2,[1,2],dimensions);
            i1o1o2X2 = myoperation(X2,[1,2,4],dimensions);
            X2 == o1X2 - i1o1X2+ i1o1o2X2;
            % objective function
            maximize(-trace(C)-4*real(trace(conj(Phi)*B*conj(H))));
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = X1;
        varargout{3} = X2;
        varargout{4} = B;
        varargout{5} = C;
    % 
    % Strategies of type iv
    % 
    case 4
        cvx_begin sdp quiet
            % variables
            variable X(dim,dim) complex semidefinite;
            variable C(components,components) complex semidefinite;
            variable B(components,dim) complex;
            expression bigmat(components+dim,components+dim);
            bigmat(1:dim,1:dim) = X;
            bigmat(dim+1:end,1:dim) = B;
            bigmat(1:dim,dim+1:end) = B';
            bigmat(dim+1:end,dim+1:end) = C;
            % restrictions
            bigmat >= 0;
            B*conj(Phi) == Phi.'*B';
            trace(X) == dimensions(2)*dimensions(4);
            i1o1o2X = myoperation(X,[1,2,4],dimensions);
            i1o1X = myoperation(X,[1,2],dimensions);
            o1i2o2X1 = myoperation(X,[2,3,4],dimensions);
            i2o2X1 = myoperation(X,[3,4],dimensions);
            o1X = myoperation(X,2,dimensions);
            o2X1 = myoperation(X,4,dimensions);
            o1o2X = myoperation(X,[2,4],dimensions);
            X == i1o1o2X- i1o1X +o1i2o2X1 -i2o2X1 + o1X+o2X1-o1o2X;
            % objective function
            maximize(-trace(C)-4*real(trace(conj(Phi)*B*conj(H))));
        cvx_end
        obj = cvx_optval;
        varargout{1} = obj;
        varargout{2} = X;
        varargout{3} = B;
        varargout{4} = C;
    otherwise
        error("Input strategy is illegal. 1 for i, 2 for ii, 3 for iii, 4 for iv.");
end
end