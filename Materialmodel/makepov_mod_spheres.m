function makepov_mod_spheres(data)

[f, p] = uiputfile('*', 'POV Filename?:');
cd(p)
fid = fopen([p f], 'w');

nsteps = size(data, 1);
spheresize = .0048*1000;
fprintf(fid, 'union { \n');

for n = 1:nsteps
    
    if n~= nsteps
        outstr = [ ' sphere{<', num2str(data(n,1)*1000), ',',...
            num2str(data(n,2)*1000), ',', num2str(data(n,3)*1000), '> ',',', num2str(spheresize),' pigment{Red}', '}'];
        fprintf(fid, '%s \n', outstr);
    else
        
        fprintf(fid, '  \n } \n');
    end
    

end


fclose(fid);

disp('Done.')