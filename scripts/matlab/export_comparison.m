function export_comparison(filename, ids, epochs, fields)
% Export explicit, identity-aligned arrays without reading a production project.
% ids: N x D numeric identities (PS coordinates or arc endpoint identities).
% epochs: numeric YYYYMMDD vector, or [] for nontemporal arrays.
% fields: struct of numeric N x T values, or N x 1 scalar metrics.
% Use a separate export for temporal and nontemporal fields.
% A matching JSON sidecar (units/reference/phase_direction) is mandatory for
% the Python comparator; it must be supplied by the reviewed production adapter.
if exist(filename, 'file')
    error('depsipy:ExistingExport', 'Choose a new export filename.');
end
if ~isnumeric(ids) || isempty(ids) || any(~isfinite(ids(:)))
    error('depsipy:InvalidIdentity', 'ids must contain finite numeric identities.');
end
if size(unique(ids, 'rows'), 1) ~= size(ids, 1)
    error('depsipy:DuplicateIdentity', 'ids must be unique.');
end
names = fieldnames(fields);
if any(ismember(names, {'ids','epochs'}))
    error('depsipy:ReservedField', 'ids and epochs are reserved field names.');
end
for k = 1:numel(names)
    a = fields.(names{k});
    if ~isnumeric(a) || size(a,1) ~= size(ids,1)
        error('depsipy:ShapeMismatch', 'Every field needs N numeric rows.');
    end
    if ~isempty(epochs) && size(a,2) ~= numel(epochs)
        error('depsipy:TimeMismatch', 'Temporal arrays must have one column per epoch.');
    end
end
fields.ids = ids;
if ~isempty(epochs)
    if ~isnumeric(epochs) || numel(unique(epochs)) ~= numel(epochs)
        error('depsipy:InvalidEpochs', 'epochs must be unique numeric dates.');
    end
    fields.epochs = epochs(:);
end
save(filename, '-struct', 'fields', '-v7');
end
