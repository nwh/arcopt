function struct2ws(s)
%struct2ws  exports all struct properties into the workspace
%
% Example:
%   s.a = 1;
%   s.b = 2;
%   struct2ws(s);
%   whos
%   % a and b are now variables in the workspace
%
% Input:
%   s = struct
%
% Result: struct properties are exported to workspace
%
%

f = fieldnames(s);
n = length(f);
for i = 1:n
  assignin('caller',f{i},s.(f{i}));
end
