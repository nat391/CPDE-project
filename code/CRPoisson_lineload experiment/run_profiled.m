function run_profiled(function_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
profile on;
function_name;
profile off;

profsave(profile('info'), 'profiler_results');

disp('Profiler-Bericht wurde in den Ordner "profiler_results" gespeichert.');
end