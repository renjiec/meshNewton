function [v, f, vt, vn] = readObjCut(filename, useCuts)
%
% obj = readObjCut(fname)
%
% This function parses wavefront object data
% It reads the mesh vertices, texture coordinates, normal coordinates
% and face definitions(grouped by number of vertices) in a .obj file 
% Optionally uses cuts from uv map
% INPUT: fname - wavefront object file full path
%      : useCuts - determine if a UV map with cuts is given and should be used
% OUTPUT: v - mesh vertices
%       : vt - texture coordinates
%       : vn - normal coordinates
%       : f - face definition assuming faces are made of of 3 vertices
%

varfacesizes = false; 

[filepath,filename,fileext] = fileparts(filename);
if isempty(fileext)
	objname = [filename '.obj'];
else
	objname = [filename fileext];
end

if isempty(filepath), filepath = '.'; end

fullfilepath = [filepath filesep objname];
% Check if obj file exists
if exist(fullfilepath, 'file') == 0
    fullfilepath = [getenv('ModelsFolder') filesep objname];
    if exist(fullfilepath, 'file') == 0
        error(['Cannot find ' objname]);
    end
end

% set up field types
str = sprintf('\n\n%s', fileread(fullfilepath));

if numel(str)>1 && str(1) == 'v'
    assert(false, 'Invalid file format!');
end

% use 'tokens' instead of 'match' to remove the text head
% a = regexp(str, '\n\s*v (.*)', 'dotexceptnewline', 'tokens');

a = regexp(str, '\n\s*v .*', 'dotexceptnewline', 'match');

% v = cell2mat(textscan(sprintf('%s\n', a{:}), 'v %f %f %f%*[^\n]'));
% v = cell2mat(textscan(strcat(a{:}), 'v %f %f %f%*[^\n]'));
v = textscan([a{:}], 'v %f %f %f%*[^\n]', 'CollectOutput', 1);
v = v{1};

a = regexp(str, '\n\s*f .*', 'dotexceptnewline', 'match');  % allow spaces before 'f ...'
if ~isempty(a)
    noslash = isempty( strfind(a{1}, '/') );
    oneslash = length(strfind(a{1}, '/'))==3;
    twoslash = length(strfind(a{1}, '/'))==6;
    if ~varfacesizes
        deg = numel( regexp(a{1}, '\d%*[^ ]*') );
        if noslash
            f = cell2mat(textscan([a{:}], ['f' repmat(' %d', 1, deg)]));
        else
            f = cell2mat(textscan([a{:}], ['f' repmat(' %d%*[^ ] ', 1, deg-1) ' %d%*[^\n]'])); % to do, take care other attributes
            if(twoslash)
              ft = cell2mat(textscan([a{:}], ['f' repmat(' %*[^/]/%d%*[^ ] ', 1, deg-1) ' %*[^/]/%d%*[^\n]']));
            elseif(oneslash)
              ft = cell2mat(textscan([a{:}], ['f %*[^/]/%d %*[^/]/%d %*[^/]/%d%*[^\n]']));
            else
              assert(false,'unsupported data');
            end
            ft = double(ft);
        end
    else
        nf = numel(a);
        f = cell(nf,1);
        for i=1:nf
            deg = numel( regexp(a{i}, '\d%*[^ ]*') );
            if noslash
                f{i} = cell2mat( textscan(a{i}, ['f' repmat(' %d', 1, deg)]) );
            else
                f{i} = cell2mat( textscan(a{i}, ['f' repmat(' %d%*[^ ] ', 1, deg-1) ' %d%*[^\n]']) ); % to do, take care other attributes
            end
        end
    end
end

vt = []; vn = []; 
a = regexp(str, '\nvn .*', 'dotexceptnewline', 'match');
if ~isempty(a)
    vn = cell2mat(textscan([a{:}], 'vn %f %f %f'));
end

a = regexp(str, '\n\s*vt .*', 'dotexceptnewline', 'match');
%a = cellfun(@(x) x(3:end), a, 'UniformOutput', false);
if ~isempty(a)
%     vt = cell2mat(textscan([a{:}], '%f %f %*[^\n]'));

%     vt = sscanf([a{:}], '%f', [numel(a),inf]);
    vt = textscan([a{:}], 'vt %f %f%*[^\n]', 'CollectOutput', 1);
    vt = vt{1};
    vt3 = zeros(size(vt,1),3);
    for i = 1:size(ft,1)
        vt3(ft(i,1),:) = v(f(i,1),:);
        vt3(ft(i,2),:) = v(f(i,2),:);
        vt3(ft(i,3),:) = v(f(i,3),:);
    end
end

if ~iscell(f)
    f = double(f);
    if min( min(f) ) == 0, f = f+1; end
else
    f = cellfun(@double, f, 'UniformOutput', false);
    if min(cellfun(@min, f)) == 0, f = cellfun(@(x) x+1, f, 'UniformOutput', false); end
end
if useCuts
  assert( size(vt,1)>1, 'No cuts in file')
  v = vt3;
  f = ft;
end
