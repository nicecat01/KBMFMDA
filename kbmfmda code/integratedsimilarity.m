function [sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km)           %分别计算miRNA相似性以及疾病相似性
sm = FS.*FSP+km.*(-(FSP-1));            %miRNA相似性
sd = SS.*SSP+kd.*(-(SSP-1));            %疾病相似性
end