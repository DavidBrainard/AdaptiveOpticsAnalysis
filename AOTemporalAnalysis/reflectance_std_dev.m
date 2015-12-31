function [ ref_stddev, ref_times ] = reflectance_std_dev( cell_times, cell_reflectance, series_length )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ref_avg = zeros( series_length,1);
ref_stddev = zeros( series_length,1);

has_content = false;
ref_times=[];
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
            ref_stddev(k) = ref_stddev(k) + ( cell_reflectance{i}(cell_times{i} == k)-ref_avg(k) ).^2;
            has_content = true;
        end
    end
    ref_stddev(k) = sqrt(ref_stddev(k)/n);
    
    if has_content
        ref_times = [ref_times; k];
    else
        ref_times = [ref_times; NaN];
    end
    
    has_content = false;
end

% ref_stddev=ref_stddev(~isnan(ref_stddev));

end

