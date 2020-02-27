[filenames, pathname] = uigetfile('*mat','Pick files','MultiSelect', 'on');


savepath = 'D:\DATA\SfN2018_ModBio\NHP-Bullet\'; saveName = 'BulletDrumTrial-';

shankA = [1:8]; shankB = [9:16]; shankC = [17:24]; shankD = [25:32];
shanks = [shankA; shankB; shankC; shankD];
probeNames = ['Probe 1'; 'Probe 2';];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];

d = [0 0.25 0.5 0.75 1.0, 1.25 1.5 1.75];

disp('Beginning Single File Parsing');
for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    disp('Plotting...')

% % 
% %     hF = figure;
% %     maximize;
    for probeNum = 1:2
        subplot(1,2,probeNum)
        for shankNum = 1:4
            tempCov = corr(filtData.shankData{probeNum}([shanks(shankNum,:)],:)'); %corr


    % zero =  diag(tempCov);
    % one =  diag(tempCov,1);
    % two =  diag(tempCov,2);
    % three =  diag(tempCov,2);
    % four = diag(tempCov, 3);
    % five = diag(tempCov, 4);
    % six = diag(tempCov, 5);
    % seven = diag(tempCov, 6);
    % eight = diag(tempCov, 7);


            fDist(shankNum, 1) = mean(diag(tempCov));
            fDistErr(shankNum,1) = std(diag(tempCov))/sqrt(length(diag(tempCov)));
            for i = 1:7
                fDist(shankNum,i+1) = mean(diag(tempCov,i));
                fDistErr(shankNum,i+1) = std(diag(tempCov,i))/sqrt(length(diag(tempCov,i)));
            end
% %             hE(probeNum, shankNum) = errorbar(d,fDist(shankNum,:), fDistErr(shankNum,:), 'LineWidth', 2);

% %             hold on;
        end

% %         ylabel('Covariance')
% %         xlabel('Distance (mm)')
% %         xticks(d)
% %         legend(shankNames);
% %         title(probeNames(probeNum,:));
        cohData(fNum).probe(probeNum).fDist = fDist;
        cohData(fNum).probe(probeNum).fDistErr = fDistErr;
        
    end
% %     suptitle([filtData.filename '- Covariance as a function of Distance: ' filtData.trialInfo]);
    clear filtData
% %     
% %     saveas(hF,['D:\Figures\SFN modBio\Figures\Filtered\CovDist\' char(erase(filenames(fNum),'.mat'))...
% %             '-filtered-CovByDist.png']);
% %     savefig(hF,['D:\Figures\SFN modBio\Figures\Filtered\CovDist\' char(erase(filenames(fNum),'.mat'))...
% %             '-filtered-CovByDist']);
     
end
% % 
% % disp('Beginning Single File Parsing');
% % for fNum = 1%1:size(filenames,2)
% %     disp('Loading file...')
% %     load(string(fullfile(pathname,filenames(fNum))));
% %     disp('Plotting...')
% % 
% % 
% %     
% %     
% %     for probeNum = 1:2
% %         subplot(1,2,probeNum)
% %         for shankNum = 1:4
% %             cD1 = filtData.shankData{probeNum}([shanks(shankNum,1)],:)';
% %             cD2 = filtData.shankData{probeNum}([shanks(shankNum,2:end)],:)';
% % 
% % 
% %            cxy = mscohere(cD1, cD2);
% %            figure;
% %            plot(cxy);
% %         end
% %     end
% % end

j = 1;
for pN = 1:2
    for sN = 1:4
        for i =1:2
            if i == 4
                continue;
            end
            tmpData(pN).shank(sN).corr(j,:) = cohData(i).probe(pN).fDist(sN,:);
            tmpError(pN).shank(sN).corr(j,:) = cohData(i).probe(pN).fDistErr(sN,:);
            j = 1+j;
        end
    end
end

hCF = figure;
for pN = 1:2
    subplot(1,2, pN)
    for sN = 1:4
        newfDist = mean(tmpData(pN).shank(sN).corr);
        newfErr = mean(tmpError(pN).shank(sN).corr);
        hE(probeNum, shankNum) = errorbar(d,newfDist, newfErr, 'LineWidth', 2);
        hold on;
    end
    ylabel('Correlation')
    xlabel('Distance (mm)')
    xticks(d)
    legend(shankNames);
    title(probeNames(pN,:));
end
suptitle(['Correlation as a function of Distance']);

saveas(hCF,['D:\Figures\SFN modBio\Figures\Filtered\CorrDist\Shaker-AggCorrByDist.png']);
savefig(hCF,['D:\Figures\SFN modBio\Figures\Filtered\CorrDist\Shaker-AggCorrByDist']);