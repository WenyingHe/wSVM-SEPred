clc
clear
cell_types={'mESC','myotube','macrophage','proB','Th','H2171','U87','MM1S'};
list=listdlg('liststring',cell_types);
cell_type=cell_types{list};%
train_model=strcat(cell_type,'_model.mat');
[FileName, PathName, FilterIndex]=uigetfile('.txt','Select file to open');
[head seq]=fastaread(FileName);
n=length(seq);
out_file=strcat(cell_type,'_predict result.fasta');
fid1=fopen(out_file,'wt');
%%
% test
Fea_3mer=Fea_kmer(FileName,3);
test_F1=mapminmax(Fea_3mer,0,1);
for k=4:7
   Fea_entro_test=['xtest_' num2str(k)];
   eval([Fea_entro_test,'=','Fea_entropy_single(Fea_kmer(FileName,k))']);
end
test_F2=[xtest_4,xtest_5,xtest_6,xtest_7];
xtest=[test_F1,test_F2];
ytest=zeros(n,1);

load(train_model);
[predict_label]=svmpredict(ytest,xtest,model,'-q');

clc;
uu=0;
for m=1:n   
    if predict_label(m)==1
        fprintf(fid1,'%s%s%s\n','>',head{1,m},'---SE');
        fprintf(fid1,'%s\n',seq{1,m});
    else
        fprintf(fid1,'%s%s%s\n','>',head{1,m},'---TE');
        fprintf(fid1,'%s\n',seq{1,m});
    end
end
