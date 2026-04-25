function run_profiled(function_name)

profile on;
feval(function_name);
profile off;

profsave(profile('info'), 'profiler_results');

disp('Profiler-Bericht wurde in den Ordner "profiler_results" gespeichert.');
end