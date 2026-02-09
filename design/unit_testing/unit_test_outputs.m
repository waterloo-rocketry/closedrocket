load('pending_tests.mat');
addpath(genpath('unit_testing'));
addpath(genpath('model'));

results = struct();

for i = 1:length(allTests)
    funcHandle = str2func(allTests(i).funcName);
    
    argList = struct2cell(allTests(i).inputs);
    
    try
        results(i).modelName = allTests(i).modelName;
        results(i).inputData = allTests(i).inputs;
        results(i).output    = funcHandle(argList{:});
    catch ME
        warning('Failed to run %s: %s', allTests(i).modelName, ME.message);
    end
end

% save results to file
save('test_results.mat', 'results');
fprintf('Results saved to test_results.mat. \n');