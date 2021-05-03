filepath = mfilename('fullpath');
filepath = strsplit(filepath,'\');
filepath{end} = 'mex';
filepath = strjoin(filepath,'\');

old_dir = cd(filepath);

c_files = dir('*.c');

for i = 1:length(c_files)
    mex(c_files(i).name);
end

cd(old_dir);
clear;