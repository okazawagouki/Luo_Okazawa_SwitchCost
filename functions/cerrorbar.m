function h = cerrorbar(X, Y, SE, varargin)
%function cerrorbar(X, Y, SE, varargin)

if length(X) ~= length(Y)
    error('length of X and Y should match');
end

if isempty(X)
    error('X is empty');
end

if size(SE,1)==1
    SE = SE';
end
if length(Y) ~= size(SE,1)
    error('length of Y and raw of SE should match');
end
if any(SE(:) < 0)
    error('SE should be all positive');
end
if size(SE,2)==1
    SE = [SE, SE];
elseif size(SE,2) > 2
    error('column of SE should be either 1 or 2');
end


h = nan(length(X),1);
hold on;
for n=1:length(X)
    h(n) = plot([1 1]* X(n), [1 1] * Y(n) + [-1 1] .* SE(n,:), varargin{:});
end



