function vout=get_augmented_toep_eig(tfcol,tfoffrow)
% this function computes the eigenvalues of the circulant matrix that
% obtained by embedding a Toeplitz matrix T into a larger circulant matrix.
%tfcol represents T(:,1), tfoffrow denotes T(1,2:end)
N=length(tfcol);
extlen=2^(ceil(log2(2*N-1)));
vout=fft([tfcol;zeros(extlen-2*N+1,1);tfoffrow(end:-1:1)]);



end