function [entro_N]=Fea_entropy_single(kmer_NN)

One_N=kmer_NN'+0.0001;
entro_N=-sum(One_N.*log2(One_N),1)';%N*1
