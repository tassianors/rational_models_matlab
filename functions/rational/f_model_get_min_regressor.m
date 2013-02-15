function min_reg = f_model_get_min_regressor(m)
%====================
%% Get max regressor id based on model information.
% m:: model structure
%
% Return max regressor id based on the model info
%====================

% Check if model parameter structure is valid
f_check_model(m);
min_reg = min(min(m.a_regress), min(m.b_regress));
end

