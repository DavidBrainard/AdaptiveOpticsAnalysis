function [ ref_variance, ref_times, ref_count ] = reflectance_pooled_variance( cell_times, cell_reflectance, series_length )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    ref_avg = zeros( series_length,1);

    has_content = false;

    ref_variance = zeros(series_length,1);
    ref_count = zeros(series_length,1);
    ref_times=[];
    
    
    
    for k=1:series_length % Plot average
        n=0;

        for i=1:length(cell_reflectance)
            if any(cell_times{i} == k) 
                ref_avg(k) = ref_avg(k) + cell_reflectance{i}(cell_times{i} == k);
                n = n+1;
            end
        end

        ref_count(k) = n;

        ref_avg(k) = ref_avg(k)/n;

        for i=1:length(cell_reflectance)
            if any(cell_times{i} == k)
                ref_variance(k) = ref_variance(k) + ( cell_reflectance{i}(cell_times{i} == k)-ref_avg(k) ).^2;
                has_content = true;
            end
        end
        
%         ref_variance(k) = ref_variance(k)/(ref_count(k)-1); % -1 For pooled variance calculation

        if has_content
            ref_times = [ref_times; k];
        else
            ref_times = [ref_times; NaN];
        end

        has_content = false;
    end

end

