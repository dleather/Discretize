function dispcalc(comm)
% function dispcalc(comm)
% disp(comm)
% evalin('caller', comm)

disp(sprintf('\n%s =\n', comm))
disp(evalin('caller', comm))