[ns_status, hFile] = ns_OpenFile('R:\data_raw\primate\NHP_DRG_Stim\Bullet_Data\Grapevine\datafile0074-01-Sorted.nev');
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);

hCF = figure;
maximize;
index = reshape(1:32, 4, 8).';

disp('Starting next Shank');
j=1;
for i = 1:32%nsFileInfo.EntityCount
    hL(i) = subplot(8,4,index(j));
    plot_spikes(hFile,i);
    disp(['Plotting ' num2str(i)]);
    j = j+1;
end
linkaxes([hL],'xy');
ylim([-500 500]);

saveas(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_Spikes.tiff']);
savefig(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_Spikes']);

pause(0.1)

close all;


hCF2 = figure;
maximize;
index2 = reshape(1:32, 4, 8).';

disp('Starting next Shank');
j=1;
for i = 33:64%nsFileInfo.EntityCount
    hL2(i) = subplot(8,4,index2(j));
    plot_spikes(hFile,i);
    disp(['Plotting ' num2str(i)]);
    j = j+1;
end
linkaxes([hL2],'xy');
ylim([-500 500]);

saveas(hCF2,['D:\Figures\SFN modBio\Figures\Filtered\Probe2_Spikes.tiff']);
savefig(hCF2,['D:\Figures\SFN modBio\Figures\Filtered\Probe2_Spikes']);

% close all
% ns_CloseFile(hFile);

% % % hCF = figure;
% % disp('Starting next Shank');
% % j = 1;
% % for i = 9:16%nsFileInfo.EntityCount
% %     hL(i) = subplot(8,8,8*(i-1)+2);
% %     plot_spikes(hFile,i);
% %     disp(['Plotting ' num2str(i)]);
% %     j=j+1;
% % end
% % linkaxes([hL],'xy');
% % ylim([-500 500]);
% % 
% % % saveas(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_ShankB.tiff']);
% % % savefig(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_ShankB']);
% % 
% % % close all
% % % ns_CloseFile(hFile);
% % 
% % % hCF = figure;
% % disp('Starting next Shank');
% % j = 1;
% % for i = 17:24%nsFileInfo.EntityCount
% %     hL(i) = subplot(8,8,8*(i-1)+3);
% %     plot_spikes(hFile,i);
% %     disp(['Plotting ' num2str(i)]);
% %     j=j+1;
% % end
% % linkaxes([hL],'xy');
% % ylim([-500 500]);
% % 
% % % saveas(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_ShankC.tiff']);
% % % savefig(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_ShankC']);
% % % 
% % % close all
% % % % ns_CloseFile(hFile);
% % % 
% % % hCF = figure;
% % disp('Starting next Shank');
% % j = 1;
% % for i = 25:32%nsFileInfo.EntityCount
% %     hL(i) = subplot(8,8,8*(i-1)+4);
% %     plot_spikes(hFile,i);
% %     disp(['Plotting ' num2str(i)]);
% %     j=j+1;
% % end
% % linkaxes([hL],'xy');
% % ylim([-500 500]);
% % 
% % saveas(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_spikes.tiff']);
% % savefig(hCF,['D:\Figures\SFN modBio\Figures\Filtered\Probe1_spikes']);
% % 
% % close all
ns_CloseFile(hFile);