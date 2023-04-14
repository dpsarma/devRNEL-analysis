%% Extract Phase from LG
figure; figure; tiledlayout(4,2); maximize;
for m = [10 2 12 4 14 6 16 8] %[2 4 6 8]
    env = envelope(abs(z(m,:)),100,'rms');
    nexttile;
    for i = [3 5 7 9 11]
        epoch = env(idx(i):idx(i+2));
        e_tim = linspace(0,100,length(epoch));
        plot(e_tim, epoch);
        hold on;
        p(m,i) = peak2peak(z(m,idx(i):idx(i+2)));
        r(m,i) = rms(z(m,idx(i):idx(i+2)));
    end
    peak(m) = max(p(m,:));
    peakrms(m) = max(r(m,:));
    title(mLabels(chan_remap(m)));
    ylabel('uV')
end
sgtitle({cell2mat(trialtype), cell2mat(filename), ['Right Heel Strike:Heel Strike']},'interpreter', 'none')

