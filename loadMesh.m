function M = loadMesh(filename)
fprintf('loading mesh %s ...', filename);
tic;
S = LoadOBJFile(filename);
M.faces = S{1}.faces'+1;
M.vertices = S{1}.vertices';
t = toc;
fprintf('finished in %f seconds.\n', t);
end