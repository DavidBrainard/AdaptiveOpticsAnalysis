function [ ref_stddev ] = reflectance_pooled_variance( cell_times, cell_reflectance, series_length )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pvar = 0;
pvar_n = 0;
for i=1:length(cell_reflectance)
    
    pvar = pvar + (length(cell_reflectance{i})-1)*var( cell_reflectance{i} );
    pvar_n = pvar_n + (length(cell_reflectance{i})-1);
end

pvar = pvar/pvar_n;
pstddev = sqrt(pvar);

% figure(10);
% hold on;
% hz=16;

ref_avg = zeros( series_length,1);
ref_stddev = zeros( series_length,1);

for k=1:series_length % Plot average
    n=0;
    
    for i=1:length(cell_reflectance)
        if any(cell_times{i} == k) 
            ref_avg(k) = ref_avg(k) + cell_reflectance{i}(cell_times{i} == k);
            n = n+1;
        end
    end
    ref_avg(k) = ref_avg(k)/n;
    
    for i=1:length(cell_reflectance)
        if any(cell_times{i} == k)
            ref_stddev(k) = ref_stddev(k) + (cell_reflectance{i}(cell_times{i} == k)-ref_avg(k)).^2;
        end
    end
    ref_stddev(k) = sqrt(ref_stddev(k)/n);
end


end

