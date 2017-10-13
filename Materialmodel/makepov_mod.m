function makepov_mod(data)

[f, p] = uiputfile('*', 'POV Filename?:');
cd(p)
fid = fopen([p f], 'w');

nsteps = size(data, 1);
spheresize = .0048*1000;
fprintf(fid, '%d \n', nsteps);

for n = 1:nsteps
    
    if n~= nsteps
        outstr = [ ' <', num2str(data(n,1)*1000), ',',...
            num2str(data(n,2)*1000), ',', num2str(data(n,3)*1000), '> ',',', num2str(spheresize),','];
        fprintf(fid, '%s \n', outstr);
    else
        outstr = [ ' <', num2str(data(n,1)*1000), ',',...
            num2str(data(n,2)*1000), ',', num2str(data(n,3)*1000), '> ',',', num2str(spheresize)];
        fprintf(fid, '%s \n', outstr);
    end
    

end


fclose(fid);

disp('Done.')