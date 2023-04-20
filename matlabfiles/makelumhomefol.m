function newf = makelumhomefol (homefol)

if ~exist(homefol)
    mkdir(homefol)
end

nom=DT4filename;
nom=nom(4:9);
newf=fullfile(homefol, nom);
tstat= mkdir(newf);
if tstat==1
cd(newf)
mkdir('matfiles')
mkdir('lumfiles')
mkdir('output')
end








end