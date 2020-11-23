figure;
imagesc(spinalmap')
colorbar
yticklabels({'L-1', 'L-2','L-3','L-4','L-5','S-1','S-2'});
xticklabels({'iliopsoas', 'adductor', 'adductor magnus', 'quadriceps', 'biceps femoris', 'gastrocnemius', 'anterior tibialis','abductor hallucis'});
xtickangle(45)
set(gca,'xaxisLocation','top')


figure;
imagesc(normalize(spinalmap','range'))
colorbar
yticklabels({'L-1', 'L-2','L-3','L-4','L-5','S-1','S-2'});
xticklabels({'iliopsoas', 'adductor', 'adductor magnus', 'quadriceps', 'biceps femoris', 'gastrocnemius', 'anterior tibialis','abductor hallucis'});
xtickangle(45)
set(gca,'xaxisLocation','top')

cmap = magma;
colormap(cmap);
figure;
Z = normalize(spinalmap','range');
Z(:, [6 7]) = Z(:, [7 6]);
imagesc(Z(:,4:7));
colorbar
yticklabels({'L-1', 'L-2','L-3','L-4','L-5','S-1','S-2'});
xticks([1 2 3 4]);
xticklabels({'quadriceps', 'biceps femoris', 'gastrocnemius', 'tibialis anterior '});
xtickangle(45)
set(gca,'xaxisLocation','top')
box off
grid off
set(gca,'TickLength',[0 0])
