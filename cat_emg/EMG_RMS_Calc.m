% EMG_RMS_Calc

%%%CONSTANTS%%%
cat={'Flahr','Gumband','Heighth','Idjit','Jims'};
type={'Fac','Mix','Inh'};
%%%%%%%%%%%%%%%

%%%INPUTS%%%
con=[6;11];
ff=3;
%aa=1;
%%%%%%%%%%%%
[responses,~,~]=xlsread('RMSValues_Final.xlsx',type{ff});
RMS_total=[];
RMS_total_week=[];
for aa=4 %cycle through cats
    %cat_response=responses(responses(:,1)==aa,:);
    weeks=[6];
    for bb=1:length(weeks) %cycle through weeks
        if weeks(bb)<10
            this_week=['0',num2str(weeks(bb))];
        else
            this_week=num2str(weeks(bb));
        end
%         cd('R:\users\ecb43\EMG Reshape\Updated Analysis Structures')
          cd('R:\users\urbinma\CatStructs')
        filename=[cat{aa},' Week ',this_week];
        %load(filename)
        for dd=1:2
            for cc=1:10
                T0=find(round(con_map(con(dd)).time*1000000)==-100000);
                T1=find(round(con_map(con(dd)).time*1000000)==0);
                T2=find(round(con_map(con(dd)).time*1000000)==20000);
                T3=find(round(con_map(con(dd)).time*1000000)==40000);
                T4=find(round(con_map(con(dd)).time*1000000)==200000);
                temp_mean=mean(con_map(con(dd)).raw(:,:,cc),2);
                RMS_0=rms(temp_mean(T0:T1));
                RMS_20=rms(temp_mean(T1:T2));
                RMS_40=rms(temp_mean(T2:T3));
                RMS_200=rms(temp_mean(T3:T4));
                RMS=horzcat(aa,weeks(bb),cc,con(dd),RMS_0,RMS_20,RMS_40,RMS_200);
                RMS_total=vertcat(RMS_total,RMS);
            end
        end
    end
end






