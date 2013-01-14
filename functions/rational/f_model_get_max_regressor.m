function max_reg = f_model_get_max_regressor(m)
%====================
%% Get max regressor id based on model information.
% m:: model structure
%
% Return max regressor id based on the model info
%====================

% Check if model parameter structure is valid
f_check_model(m);
max_reg = max(max(abs(m.a_regress)), max(abs(m.b_regress)));
