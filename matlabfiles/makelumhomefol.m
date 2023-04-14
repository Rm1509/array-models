function newf = makelumhomefol (homefol)

cd(homefol)
nom=DT4filename;
nom=nom(4:9);
newf=strcat(homefol, '/',nom);
tstat= mkdir(newf);
if tstat==1
cd(newf)
mkdir('matfiles')
mkdir('lumfiles')
mkdir('output')
end








end