%%
clear;
clc;
close all;
load('fp_database.dat','-mat');

for kk=1:fp_number
%     s1=0;
%     s2=0;
%     s3=0;
%     s4=0;
    namefile=char(namefile_vector(kk));
    S=regexp(namefile,'.jpg','split');
    str=char(S(1));
    outfeature=['f' str '.txt'];
    fpathname='H:\LAB\Dataset\';
    fd=fopen(strcat(fpathname,outfeature),'wt');
    for n = 1:8
        feature=data{kk,1}{1,n}(:,1)';
        for k = 1:64
            if feature(k) >= 300
                a='1';
%                 s1=s1+1;
            elseif feature(k)>=200
                a='2';
%                 s2=s2+1;
            elseif feature(k)>=100
                a='3';
%                 s3=s3+1;
%             elseif feature(k)>=50
%                 a='4';
%                 s4=s4+1;
            else
                a='5';
%                 s4=s4+1;
            end
            fprintf(fd,'%s',a);
        end
    end  
    fclose(fd);
%     disp(s1);
%     disp(s2);
%     disp(s3);
%     disp(s4);
%     disp('E');
end