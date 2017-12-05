function [ ref_variance, ref_times, ref_count ] = reflectance_pooled_variance( cell_times, cell_reflectance, series_length )
% [ ref_variance, ref_times, ref_count ] = reflectance_pooled_variance( cell_times, cell_reflectance, series_length )
%   This function calculates the pooled variance of a group of signals, where each time point
%   may have a different number of signals contributing to it.
%
% @params:
%    cell_times: A 1xN cell array of frame indexes. Each cell contains to the frame times 
%                from a single photoreceptor signal.
%
%    cell_reflectance: A 1xN cell array of photoreceptor reflectances.
%                      Each cell contains the reflectance signal from a single photoreceptor.
%
%    series_length: The maximum length that the cell_reflectance or cell_times can be.
%
% @outputs:
%
%    ref_variance: An array containing the variance at each time point that we have data.
%
%    ref_times: An array containg the frame indexes of the ref_variance
%
%    ref_count: An array as long as ref_times where each array entry contains the number of signals that
%               contributed to each time point.

    ref_avg = zeros( series_length,1);

    has_content = false;

    ref_variance = zeros(series_length,1);
    ref_count = zeros(series_length,1);
    ref_times=[];
    
    
    
    for k=1:series_length % Plot average
        n=0;

        for i=1:length(cell_reflectance)
            if any(cell_times{i} == k) &&...
           ~isinf(cell_reflectance{i}(cell_times{i} == k)) &&...
           ~isnan(cell_reflectance{i}(cell_times{i} == k)) 
                ref_avg(k) = ref_avg(k) + cell_reflectance{i}(cell_times{i} == k);
                n = n+1;
            end
        end

        ref_count(k) = n;

        ref_avg(k) = ref_avg(k)/n;

        for i=1:length(cell_reflectance)
            if any(cell_times{i} == k) &&...
           ~isinf(cell_reflectance{i}(cell_times{i} == k)) &&...
           ~isnan(cell_reflectance{i}(cell_times{i} == k))
                ref_variance(k) = ref_variance(k) + ( cell_reflectance{i}(cell_times{i} == k)-ref_avg(k) ).^2;
                has_content = true;
            end
        end
        
%         ref_variance(k) = ref_variance(k)/(ref_count(k)-1); % -1 For pooled std dev calculation

        if has_content
            ref_times = [ref_times; k];
        else
            ref_times = [ref_times; NaN];
        end

        has_content = false;
    end

end

