function ptx_edgeVectorization(ptx_adjmatpath, SC_outpath, volNormSC_outpath, connProb_outpath, meanLength_outpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in probabilistic connectivity matrix %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptx_adjmatpath
load(ptx_adjmatpath);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Streamline Count %%%
%%%%%%%%%%%%%%%%%%%%%%%%
SC_connVec=squareform(A_sc_und);
dlmwrite(SC_outpath, SC_connVec,'-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Volume-normalized Streamline count %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volNormSC_connVec=squareform(A_volNorm_sc_und);
dlmwrite(volNormSC_outpath, volNormSC_connVec,'-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Connectivity Probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connProb_connVec=squareform(A_volNorm_sc_und);
dlmwrite(connProb_outpath, connProb_connVec,'-append')

%%%%%%%%%%%%%%%%%%%%
%%% Fiber Length %%%
%%%%%%%%%%%%%%%%%%%%
meanLength_connVec=squareform(A_length_und);
dlmwrite(meanLength_outpath, meanLength_connVec,'-append')

end
