clc;clear;
str = load('.\1.miRNA-disease association data\knowndiseasemirnainteraction.txt');
[~,disease]=xlsread('.\1.miRNA-disease association data\disease number.xlsx');
[~,miRNA]=xlsread('.\1.miRNA-disease association data\miRNA number.xlsx');
% nd:the number of diseases
% nm:the number of miRNAs
% pp:the number of known diseae-miRNA associations
nd = max(str(:,2));
nm = max(str(:,1));
[pp,~] = size(str);
% FS:the function similarity between m(i) and m(j)
% FSP:miRNA functional similarity weighting matrix
% SS:the semantic similarity between d(i) and d(j).
% SSP:Disease semantic similarity weighting matrix 1
%
FS = load('.\4.miNA functional similarity\MiRNA functional similarity.txt');
FSP = load('.\4.miNA functional similarity\miRNA functional similarity weighting matrix.txt');         
SS1 = load('.\2.disease semantic similarity 1\Disease semantic similarity model 1.txt');
SS2 = load('.\3.disease semantic similarity 2\Disease semantic similarity model 2.txt');
SS = (SS1+SS2)/2;
SSP = load('.\2.disease semantic similarity 1\Disease semantic similarity weighting matrix 1.txt');
%interaction: adajency matrix for the disease-miRNA association network
%interaction(i,j)=1 means miRNA j is related to disease i
interaction = zeros(nd,nm);
for i = 1:pp
    interaction(str(i,2),str(i,1)) = 1;
end
[kd,km] = gaussiansimilarity(interaction,nd,nm);                   %calculating Gaussian similarity
[sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km);               %Calculate miRNA similarities and disease similarities

rand('state', 1606); %#ok<RAND>
randn('state', 1606); %#ok<RAND>
Kx=sd;
Kz=sm;
Y=interaction;
state = kbmf_regression_train(Kx, Kz, Y, 15);
prediction = kbmf_regression_test(Kx, Kz, state);
ka=state.Ax.mu';kz=state.Az.mu';
[score]=(ka*sd)' *(kz*sm);                                   
%Use the KBMF algorithm to calculate the score
%everyresult(disease,miRNA,interaction,score);             %Calculate and rank the predicted scores for each disease and write the result to excel
allresult(disease,miRNA,interaction,score);               %Calculate all the potential association scores and write the result to excel



