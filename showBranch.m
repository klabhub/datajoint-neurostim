function showBranch(repoPath)

oldDir = pwd;
cd(repoPath);
% Run git command in that folder

[status, branch] = system('git rev-parse --abbrev-ref HEAD');
cd(oldDir);

if status == 0
    branch = strtrim(branch); % remove whitespace/newline
    fprintf('Repository at %s is on branch: %s\n', repoPath, branch);
else
    warning('Could not determine Git branch for %s', repoPath);
end
end