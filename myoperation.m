function res_X = myoperation(X,sys,dimensions)
% MYOPERATION  function for operation O_X = I/d_O\otimes partialtrace(X)
% sys is the dimension to be operated
% dimensions is the dimension of X
len = length(sys);
res_X = X;
for i = 1:len
    idx = sys(i);
    res_X = PartialTrace(res_X,idx,dimensions);
    if idx == 1
        res_X = kron(eye(dimensions(idx))/dimensions(idx),res_X);
    elseif idx == length(dimensions)
        res_X = kron(res_X,eye(dimensions(idx))/dimensions(idx));
    else
        res_X = kron(eye(dimensions(idx))/dimensions(idx),res_X);
        res_X = Swap(res_X,[1,2],[dimensions(idx),prod(dimensions(1:idx-1)),dimensions(idx+1:end)]);
    end
end
end