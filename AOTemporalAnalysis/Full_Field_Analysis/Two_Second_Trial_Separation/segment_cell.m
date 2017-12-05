function [ pad_roi, adj_mat, rowcol_inds ] = segment_splitcell( roi )
% [ pad_roi, adj_mat, rowcol_inds ] = segment_splitcell( roi )
%   This function creates the cost matrix for a given roi.


% Zero pad the top and bottom
roi = [zeros(1,size(roi,2)); roi; zeros(1,size(roi,2))];
% 255-pad the top and bottom
% roi = [max(roi(:))*ones(1,size(roi,2)); roi; max(roi(:))*ones(1,size(roi,2))];


adj_mat = single(zeros(length(roi(:)),length(roi(:))));

[indi, indj] = ind2sub(size(roi), 1:length(roi(:)));

rowcol_inds = [indi' indj'];

mincost = 0.00001;

allneighbors = zeros(5,1);

for i=1:length(roi(:))
       
    IoI_value = roi(i);
    
    % 5-connected neighbors (down, diagonal, left-right only)
    allneighbors = zeros(5,1);
    left  = [];
    right = [];
    
    % Left (if possible)
%     if (rowcol_inds(i,2) > 1)
%         
%         adj_mat(i,i-size(roi,1)) = -IoI_value + roi(i-size(roi,1));
%         
%         if adj_mat(i,i-size(roi,1)) == 0 || rowcol_inds(i,1) == 1
%             adj_mat(i,i-size(roi,1)) = mincost;
%         end
%         
%         allneighbors(1) = adj_mat(i,i-size(roi,1));
%         
%     end
%     
%     % Right (if possible)
%     if (rowcol_inds(i,2) < size(roi,2))
%         
%         adj_mat(i,i+size(roi,1)) = -IoI_value + roi(i+size(roi,1));
%         
%         if adj_mat(i,i+size(roi,1)) == 0 || rowcol_inds(i,1) == 1
%            adj_mat(i,i+size(roi,1)) = mincost;
%         end
%         
% %         allneighbors = [allneighbors adj_mat(i,i+size(roi,1))];
%         allneighbors(2) = adj_mat(i,i+size(roi,1));
%     end
    
    % Diagonal-right (if possible)
    if (rowcol_inds(i,2) < size(roi,2)) && (rowcol_inds(i,1) < size(roi,1))
        
        adj_mat(i,i+size(roi,1)+1) = -IoI_value + roi(i+size(roi,1)+1);
        
        if adj_mat(i,i+size(roi,1)+1) == 0 || rowcol_inds(i,1) == 1
            adj_mat(i,i+size(roi,1)+1) = mincost;
        end
        
        right = adj_mat(i,i+size(roi,1)+1);
%         allneighbors = [allneighbors adj_mat(i,i+size(roi,1)+1)];
        allneighbors(3) = adj_mat(i,i+size(roi,1)+1);
    end
    
    % Diagonal-left (if possible)
    if (rowcol_inds(i,2) >1) && (rowcol_inds(i,1) < size(roi,1))
        
        adj_mat(i,i-size(roi,1)+1) = -IoI_value + roi(i-size(roi,1)+1);
        
        if adj_mat(i,i-size(roi,1)+1) == 0 || rowcol_inds(i,1) == 1
            adj_mat(i,i-size(roi,1)+1) = mincost;
        end
        left = adj_mat(i,i-size(roi,1)+1);
%         allneighbors = [allneighbors adj_mat(i,i-size(roi,1)+1)];
        allneighbors(4) = adj_mat(i,i-size(roi,1)+1);
    end
    
    % Down
    if (rowcol_inds(i,1) < size(roi,1))
        
        % If there is a large difference between the two directions
        % (left/right), then encourage the path to keep going down the middle of them.
        % In this case we want brighter lines to the left of the curve.
        if ~isempty(left) && ~isempty(right)
           adj_mat(i,i+1) = (-IoI_value + roi(i+1)) - left-right*.5;
        else
            if isempty(left)
                left = 0;
            end
            if isempty(right)
                right = 0;
            end
            % If we're on an edge, get away from it!
            adj_mat(i,i+1) = (-IoI_value + roi(i+1)); % + abs(left+right)*2;
        end
        
        if adj_mat(i,i+1) == 0 || rowcol_inds(i,1) == 1
           adj_mat(i,i+1) = mincost;
        end
        
%         allneighbors = [allneighbors adj_mat(i,i+1)];
        allneighbors(5) = adj_mat(i,i+1);
    end
        
    % If there are any negative weights, add the minimum
    minofeach = min(allneighbors(allneighbors~=0));
    
    %Down
    if (rowcol_inds(i,1) < size(roi,1)) && adj_mat(i,i+1) ~= mincost;
        adj_mat(i,i+1) = adj_mat(i,i+1) - minofeach + mincost;
    end
%     %Left
%     if (rowcol_inds(i,2) > 1) && adj_mat(i,i-size(roi,1)) ~= mincost;
%         adj_mat(i,i-size(roi,1)) = adj_mat(i,i-size(roi,1)) - minofeach + mincost;
%     end
%     %Right
%     if (rowcol_inds(i,2) < size(roi,2)) && adj_mat(i,i+size(roi,1)) ~= mincost;
%         adj_mat(i,i+size(roi,1)) = adj_mat(i,i+size(roi,1)) - minofeach + mincost;
%     end
    % Diagonal-right
    if (rowcol_inds(i,2) < size(roi,2)) && (rowcol_inds(i,1) < size(roi,1)) && adj_mat(i,i+size(roi,1)+1) ~= mincost;
        adj_mat(i,i+size(roi,1)+1) = adj_mat(i,i+size(roi,1)+1) - minofeach + mincost;
    end
    % Diagonal-left
    if (rowcol_inds(i,2) >1) && (rowcol_inds(i,1) < size(roi,1)) && adj_mat(i,i-size(roi,1)+1) ~= mincost;
        adj_mat(i,i-size(roi,1)+1) = adj_mat(i,i-size(roi,1)+1) - minofeach + mincost;
    end
end
%     adj_mat = abs(adj_mat);
    pad_roi = roi;

end

