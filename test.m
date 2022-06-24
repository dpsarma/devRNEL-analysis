for e = size(emg,1)
    asd = emg(e,:);
    emgSteps.(e).steps = cell2mat(arrayfun(@(x) asd(starts_L(x):starts_L(x+1)), 1:length(starts_L)-1, 'UniformOutput',false)');
    clear asd
end