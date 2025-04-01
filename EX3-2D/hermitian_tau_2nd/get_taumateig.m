function vout=get_taumateig(tfcol)
vout=tfcol-[tfcol(3:end);0;0];
vout=dst(vout);
m=length(tfcol);
theta=pi/(m+1);
smatfcol=sin((1:m)*theta).';
vout=vout./smatfcol;

end