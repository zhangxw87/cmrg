% Repeat Example 3 100 times
clear; clc;

if matlabpool('size') == 0
    matlabpool open
end

parfor loopvar = 1:100
    Example3(loopvar);
end

matlabpool close

exit
