function f_check_model(model)
%% Model
% m.dim:: Total dimension (num+den dimensions)
% m.n_dim:: dimension of numerator
% m.a_exp::exponential coef from the model
% m.a_signal_type:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m.a_exp.
% m.a_regress:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%
reqFields={'n_dim' 'dim' 'a_exp' 'a_signal_type' 'a_regress' 'b_exp' 'b_signal_type' 'b_regress'};
    
for i = 1:numel(reqFields) 
    if isfield(model, (reqFields{i})) == 0
        error('Model did not have #%s field', reqFields{i});
    end
end

if model.n_dim > model.dim
    error('Invalid model Dimension #%s must be smaller than #%s', reqFields{1}, reqFields{2});
end

if max(size(model.a_exp)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{3}, reqFields{2});
end

if max(size(model.a_signal_type)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{4}, reqFields{2});
end

if max(size(model.a_regress)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{5}, reqFields{2});
end

if max(size(model.b_exp)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{6}, reqFields{2});
end

if max(size(model.b_signal_type)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{7}, reqFields{2});
end

if max(size(model.b_regress)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{8}, reqFields{2});
end

end