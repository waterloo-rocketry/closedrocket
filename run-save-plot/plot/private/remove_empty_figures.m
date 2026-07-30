function figs = remove_empty_figures(figs)
%REMOVE_EMPTY_FIGURES Drop fields for dashboards that had no signals.

    names = fieldnames(figs);
    for i = 1:numel(names)
        if isempty(figs.(names{i}))
            figs = rmfield(figs, names{i});
        end
    end
end
