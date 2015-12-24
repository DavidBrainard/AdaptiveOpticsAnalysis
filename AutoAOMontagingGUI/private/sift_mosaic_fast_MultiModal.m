function [bestH, numOkMatches_all, numMatches_all]= sift_mosaic_fast_MultiModal(im1, im2, saveMatchesName,saveFlag,f1,d1,f2,d2,TransType)
% SIFT_MOSAIC Demonstrates matching two images using SIFT and RANSAC
%
%   SIFT_MOSAIC demonstrates matching two images based on SIFT
%   features and RANSAC and computing their mosaic.
%
%   SIFT_MOSAIC by itself runs the algorithm on two standard test
%   images. Use SIFT_MOSAIC(IM1,IM2) to compute the mosaic of two
%   custom images IM1 and IM2.

% AUTORIGHTS

% --------------------------------------------------------------------
%                                                         SIFT matches
% --------------------------------------------------------------------
MN = size(f1,1);

X1 = [];
X2 = [];
matches = cell(MN,1);
numMatches = zeros(MN,1);
for m = 1:MN
    [matches_m, scores] = vl_ubcmatch(d1{m},d2{m}) ;
    X1_m = f1{m}(1:2,matches_m(1,:)) ;
    X2_m = f2{m}(1:2,matches_m(2,:)) ;
    %check for duplicate matches
    %check X1 first
    [uX1_m, IA, IC1] = unique(round(X1_m'),'rows');
    X1_m = X1_m(1:2,IA');
    X2_m = X2_m(1:2,IA');
    matches_m = matches_m(:,IA');
    
    %check X2 next
    [uX2_m, IA, IC1] = unique(round(X2_m'),'rows');
    X1_m = X1_m(1:2,IA');
    X2_m = X2_m(1:2,IA');
    matches_m = matches_m(:,IA');
    
    matches{m} = matches_m;
    numMatches(m) = size(X1_m,2);
    X1=[X1 X1_m];
    X2=[X2 X2_m];
end


X1(3,:) = 1 ;
X2(3,:) = 1 ;


%numMatches = size(matches,2) ;
numMatches_all = size(X1,2);

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------

clear H score ok ;
bestScore = -1;
bestH = zeros(3,3);
bestOK_all = zeros(size(X1));
for t = 1:5000
  % estimate model
  subset = vl_colsubset(1:numMatches_all, 3) ;
    
  H = eye(3,3);
  if (TransType == 0)
      %Score Using Only Translation
      C1 = mean(X1(1:2,subset),2);
      C2 = mean(X2(1:2,subset),2);
      trans = C2 - C1;
      H(1,3) = trans(1);
      H(2,3) = trans(2);

  elseif (TransType == 1)
      %Score Using Only rotation + translation
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false); 
      H = [trans.T' mean(trans.c,1)'; 0 0 1];  

  elseif (TransType == 2)
      %Score Using rotation + translation + scalling
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', true, 'reflection', false); 
      H = [trans.b*trans.T' mean(trans.c)'; 0 0 1];  
  
  elseif(TransType == 3)
      %Score using homography
       A = [] ;
       for i = subset
          A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
       end
       [U,S,V] = svd(A) ;
       H = reshape(V(:,9),3,3) ;
  elseif(TransType == 4)
       %individual scale
      [dist,Z,trans]  = procrustes(X2(1:2,subset)',X1(1:2,subset)','scaling', false, 'reflection', false); 
      H = [trans.T' mean(trans.c)'; 0 0 1];
      a=X2(1,subset)*pinv(Z(:,1)');
      b=X2(2,subset)*pinv(Z(:,2)');
      H = [a 0 0; 0 b 0; 0 0 1] * H;
  end
   
   
  X2_ = H * X1 ;
  du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
  dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
  ok = (du.*du + dv.*dv) < 6*6 ;
  score = sum(ok) ;

  if(TransType == 4)
      a=X2(1,ok)*pinv(X2_(1,ok));
      b=X2(2,ok)*pinv(X2_(2,ok));
      H = [a 0 0; 0 b 0; 0 0 1] * H;

      X2_ = H * X1 ;
      du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
      dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
      ok = (du.*du + dv.*dv) < 6*6 ;
      score = sum(ok) ;
  end
  
  if(score > bestScore)
    bestScore = score;
    bestH = H;
    bestOK_all = ok;  
  end

end

numOkMatches_all = sum(bestOK_all);
numOkMatches = zeros(MN,1);
%separate valid matches by modality (visualization purpose only)
bestOK = cell(MN,1);
offset = 0;
for m = 1:MN
   bestOK{m} = bestOK_all(offset+1:offset+numMatches(m));
   numOkMatches(m) = sum(bestOK{m}); 
   offset = offset + numMatches(m);
end


% --------------------------------------------------------------------
%                                                  Optional refinement
% --------------------------------------------------------------------

% function err = residual(H)
%  u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
%  v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
%  d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
%  du = X2(1,ok) - u ./ d ;
%  dv = X2(2,ok) - v ./ d ;
%  err = sum(du.*du + dv.*dv) ;
% end
% 
% if exist('fminsearch') == 2
%   H = H / H(3,3) ;
%   opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
%   H(1:8) = fminsearch(@residual, H(1:8)', opts) ;
% else
%   warning('Refinement disabled as fminsearch was not found.') ;
% end

% --------------------------------------------------------------------
%                                                         Show matches
% --------------------------------------------------------------------

% figID=figure(1) ; clf ;
%title(sprintf('%d total tentative matches', numMatches_all)) ;
%title(sprintf('%d (%.2f%%) total inliner matches out of %d', ...
%              sum(bestOK_all), ...
%              100*sum(bestOK_all)/numMatches_all, ...
%              numMatches_all)) ;


for m = 1:MN
dh1 = max(size(im2{m},1)-size(im1{m},1),0) ;
dh2 = max(size(im1{m},1)-size(im2{m},1),0) ;

% subplot(2,MN,m) ;
% imagesc([padarray(im1{m},dh1,'post') padarray(im2{m},dh2,'post')]) ;
% o = size(im1{m},2) ;
% line([f1{m}(1,matches{m}(1,:));f2{m}(1,matches{m}(2,:))+o], ...
%      [f1{m}(2,matches{m}(1,:));f2{m}(2,matches{m}(2,:))]) ;
% title(sprintf('%d tentative matches', numMatches(m))) ;
% axis image off ;
% 
% subplot(2,MN,MN+m) ;
% imagesc([padarray(im1{m},dh1,'post') padarray(im2{m},dh2,'post')]) ;
% o = size(im1{m},2) ;
% line([f1{m}(1,matches{m}(1,bestOK{m}));f2{m}(1,matches{m}(2,bestOK{m}))+o], ...
%      [f1{m}(2,matches{m}(1,bestOK{m}));f2{m}(2,matches{m}(2,bestOK{m}))]) ;
% title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
%               sum(bestOK{m}), ...
%               100*sum(bestOK{m})/numMatches(m), ...
%               numMatches(m))) ;
% axis image off ;
end

% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
% 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% totaltext = sprintf('%d (%.2f%%) total inliner matches out of %d in all modalities', ...
%               sum(bestOK_all), ...
%               100*sum(bestOK_all)/numMatches_all, ...
%               numMatches_all) ;
% 
% text(0.5, 1,totaltext,'HorizontalAlignment' ...
% ,'center','VerticalAlignment', 'top')
% 
% 
% 
% 
% drawnow ;
% colormap('gray')
% if(saveFlag)
%     saveas(figID,saveMatchesName)
% end

% --------------------------------------------------------------------
%                                                               Mosaic
% --------------------------------------------------------------------
% 

for m = 1:MN
box2 = [1  size(im2{m},2) size(im2{m},2)  1 ;
        1  1           size(im2{m},1)  size(im2{m},1) ;
        1  1           1            1 ] ;
box2_ = inv(H) * box2 ;
box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
ur = min([1 box2_(1,:)]):max([size(im1{m},2) box2_(1,:)]) ;
vr = min([1 box2_(2,:)]):max([size(im1{m},1) box2_(2,:)]) ;

[u,v] = meshgrid(ur,vr) ;
im1_ = vl_imwbackward(im2double(im1{m}),u,v) ;

z_ = H(3,1) * u + H(3,2) * v + H(3,3);
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
im2_ = vl_imwbackward(im2double(im2{m}),u_,v_) ;

mass = ~isnan(im1_) + ~isnan(im2_) ;
im1_(isnan(im1_)) = 0 ;
im2_(isnan(im2_)) = 0 ;
mosaic = (im1_ + im2_) ./ mass ;

% figure(m+1) ; clf ;
% imagesc(mosaic) ; axis image off ;
% colormap('gray')
% title('Mosaic') ;
end
