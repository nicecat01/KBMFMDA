% Mehmet Gonen (mehmet.gonen@gmail.com)

function prediction = kbmf_regression_test(Kx, Kz, state)
    directory = fileparts(mfilename('fullpath'));
    addpath([directory, '/kbmf1k1k']);
    prediction = state.parameters.test_function(Kx, Kz, state);
end
