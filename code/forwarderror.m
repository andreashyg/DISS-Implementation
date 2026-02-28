function [forwErr] = forwarderror(x_1,x_2)
%FORWARDERROR norm of difference divided by norm of x_2
arguments (Input)
    x_1
    x_2
end

arguments (Output)
    forwErr
end

forwErr = norm(x_1 - x_2, 1) / norm(x_2, 1);
end