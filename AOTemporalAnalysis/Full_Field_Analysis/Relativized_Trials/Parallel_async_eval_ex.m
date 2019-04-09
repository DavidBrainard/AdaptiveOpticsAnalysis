

p = gcp();
% To request multiple evaluations, use a loop.
for idx = 1:10
  f(idx) = parfeval(p,@magic,1,idx); % Square size determined by idx
end
% Collect the results as they become available.
magicResults = cell(1,10);
for idx = 1:10
  % fetchNext blocks until next results are available.
  [completedIdx,value] = fetchNext(f);
  magicResults{completedIdx} = value;
  fprintf('Got result with index: %d.\n', completedIdx);
end