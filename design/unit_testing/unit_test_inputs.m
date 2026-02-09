%% Random Input Generation Function
function testSet = generateInputs(modelConfig)
    testSet = struct();

    for p = 1:length(modelConfig.params)
        currParam = modelConfig.params{p};

        if isfield(currParam, 'value')
            val = currParam.value;
        else
            val = currParam.range(1) + (currParam.range(2) - currParam.range(1)) .* rand(currParam.size);
        end

        testSet.(currParam.name) = val;
    end
end

% models to run unit tests for
models(1).name = 'AccelModel';
models(1).funcName = 'model_acceleration';
models(1).params = {
    struct('name', 'x', 'size', [10,1], 'range', [-10,10])
    struct('name', 'IMU_1', 'size', [3,1], 'range', [-10,10])
    struct('name', 'IMU_2', 'size', [3,1], 'range', [-10,10])
    struct('name', 'sensor_select', 'value', [1,1])
};

models(2).name = 'AltitudeDataModel';
models(2).funcName = 'model_altdata';
models(2).params = {
     struct('name', 'pressure', 'size', [1,1], 'range', [0,10])
};

models(3).name = 'EncoderMeasModel';
models(3).funcName = 'model_meas_encoder';
models(3).params = {
    struct('name', 't', 'size', [1,1], 'range', [-10,10])
    struct('name', 'x', 'size', [13,1], 'range', [-10,10])
    struct('name', 'bias', 'size', [1,1], 'range', [-10,10])
};

allTests = struct('modelName', "", 'funcName', "", 'inputs', []);

for i = 1:length(models)
    testInputs = generateInputs(models(i)); 
    
    allTests(i).modelName = models(i).name;
    allTests(i).funcName  = models(i).funcName;
    allTests(i).inputs    = testInputs;
end

% save to file so output generator can use inputs later
save('pending_tests.mat', 'allTests');
fprintf('Saved %d test cases to pending_tests.mat\n', length(allTests));