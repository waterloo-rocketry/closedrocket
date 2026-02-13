%% Random Input Generation Function
function testSet = generateInputs(modelConfig)
    testSet = struct();

    for p = 1:length(modelConfig.params)
        currParam = modelConfig.params{p};

        if isfield(currParam, 'value')
            val = currParam.value;
        elseif isfield(currParam, 'intervals')
            val = zeros(currParam.size);

            for r = 1:size(currParam.range, 1)
                index = currParam.intervals{r};
                rangeMin = currParam.range(r, 1);
                rangeMax = currParam.range(r, 2);
                val(index) = rangeMin + (rangeMax - rangeMin) .* rand(size(val(index)));
            end
        else
            val = currParam.range(1) + (currParam.range(2) - currParam.range(1)) .* rand(currParam.size);
        end

        testSet.(currParam.name) = val;
    end
end

% models to run unit tests for
models(1).funcName = 'model_acceleration';
models(1).params = {
    struct('name', 'x', 'size', [10,1], 'range', [-10,10])
    struct('name', 'IMU_1', 'size', [3,1], 'range', [-10,10])
    struct('name', 'IMU_2', 'size', [3,1], 'range', [-10,10])
    struct('name', 'sensor_select', 'value', [1,1])
};

models(2).funcName = 'model_altdata';
models(2).params = {
     struct('name', 'pressure', 'size', [1,1], 'range', [0,10])
};

models(3).funcName = 'model_meas_encoder';
models(3).params = {
    struct('name', 't', 'size', [1,1], 'range', [-10,10])
    struct('name', 'x', 'size', [13,1], ...
           'range', [-10,10; -5,5], ...
           'intervals', {{1:4, 5:13}})
    struct('name', 'bias', 'size', [1,1], 'range', [-10,10])
};

% number of batches of inputs
batchSize = 3;
allTests = struct('funcName', "", 'inputs', []);
testsIndex = 1;

for i = 1:length(models)
    for j = 1:batchSize
        testInputs = generateInputs(models(i));

        allTests(testsIndex).funcName = models(i).funcName;
        allTests(testsIndex).inputs = testInputs;

        testsIndex = testsIndex + 1;
    end
end

% save to file so output generator can use inputs later
save('pending_tests.mat', 'allTests');
fprintf('Saved %d test cases to pending_tests.mat\n', length(allTests));