function [ r ] = control_fit( params, data, timeBase )

%Exp
% prediction = params(1) + params(2)*params(3).^( timeBase*params(4) );

%Logistic
prediction = params(1) + ( params(2) ./ (1+exp(-params(3).*(timeBase-params(4)))) );


diff = (data-prediction).^2;
r = 100*sqrt(sum(diff(~isnan(diff)))/length(diff(~isnan(diff))));

end

