function [files] = get_files(head_dir,model,var,realm)

if ~strcmp(head_dir(end),'/')
    head_dir = [head_dir,'/'];
end

if nargin>3
    if strcmp(realm,'none') % CMIP6
        full_dir = [head_dir,model,'/',var,'/'];
    else
        full_dir = [head_dir,realm,'/',var,'/',model,'/'];
    end
elseif nargin > 2 % CMIP6
    full_dir = [head_dir,model,'/',var,'/'];
elseif nargin>1
    full_dir = [head_dir,model,'/'];
else
    full_dir = head_dir;
end

if ~exist(full_dir, 'dir')
    error(['Error: Directory ',full_dir,' does not exist.'])
end
info = dir(full_dir);

istart = 1;
while length(info(istart).name) < 3
    istart = istart + 1;
end

while strfind(info(istart).name,'DS_Store')
    istart = istart + 1;
end

if strcmp(info(istart).name,'r1i1p1')
    full_dir = [full_dir,'r1i1p1/'];
    info = dir(full_dir);
    istart = 1;
    while length(info(istart).name) < 3
        istart = istart + 1;
    end
    while strfind(info(istart).name,'DS_Store')
        istart = istart + 1;
    end
end

for i = istart:length(info)
    j = i - istart + 1;
    files{j} = [full_dir,info(i).name];
end
