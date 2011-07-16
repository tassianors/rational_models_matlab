function f_check_model(model)

mFields = fieldnames(model);
reqFields={'n_dim' 'dim' 'texp' 'yu' 'regr'};
    
% if numel(mFields) < numel(reqFields)
%     error('Model must have at least %d fields', numel(reqFields));
% end

for i = 1:numel(reqFields) 
    if isfield(model, (reqFields{i})) == 0
        error('Model did not have #%s field', reqFields{i});
    end
end

if model.n_dim > model.dim
    error('Invalid model Dimension #%s must be smaller than #%s', reqFields{1}, reqFields{2});
end

if max(size(model.texp)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{3}, reqFields{2});
end

if max(size(model.yu)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{4}, reqFields{2});
end

if max(size(model.regr)) ~= model.dim
    error('Field #%s size must be the same as #%s', reqFields{5}, reqFields{2});
end

end