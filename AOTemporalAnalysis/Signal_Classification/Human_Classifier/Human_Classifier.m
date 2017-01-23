function []=Human_Classifier()
clear;
close all force;

basePath = which('Human_Classifier.m');
basePath = basePath(1:end-(length('Human_Classifier.m')+1));
load(fullfile(basePath,'data.mat'), '-mat','x04','x45','x92');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Contents far below. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Do not scroll down, just press go.  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










































































































































data = [x04; x45; x92];
order = randperm(length(data));

ranking = [repmat('s',25,1); repmat('n',50,1)];
humanranking = repmat('',length(data),1);

fig = figure(1);
plotax = plot(1,1); hold on;
plot(66:99,5.75,'r*'); hold off;
set(fig,'KeyPressFcn',@keypress);
xlabel('Frame #'); ylabel('Standardized Reflectance')
title('Press ''S'' for a profile that shows a response, or ''N'' for a profile that shows no response');
axis([0 250 -6 6]);
for i=1:length(order)
   
    set(plotax ,'XData',data{order(i),1});
    set(plotax ,'YData',data{order(i),2});
    
    uiwait;
end

agreement = zeros(length(ranking),1);
for j=1:length(ranking)
    agreement(j) = strcmpi(ranking(j),humanranking(j));
end

disp(['You agreed with ' num2str( 100*sum(agreement)./length(agreement) ) '% of the classifications made using the wavelet/thresholding method.']);
save(fullfile(basePath,'agreement.mat'),'agreement','order');
function keypress(src, evt)

    if strcmpi(evt.Key,'s') || strcmpi(evt.Key,'n')
        humanranking(order(i)) = evt.Key;
        uiresume;
    end
    
end

end

