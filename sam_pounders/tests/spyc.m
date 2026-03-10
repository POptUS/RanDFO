function spyc(sA,cmap,pb)
%SPYC Visualize sparsity pattern with color-coded scale.
%   SPYC(S) plots the color-coded sparsity pattern of the matrix S.
%
%   SPYC(S,CMAP) plots the sparsity pattern of the matrix S USING 
%                    COLORMAP CMAP.
%
%   SPYC(S,CMAP,PB) allows turning off the display of a colorbar by passing
%                   flag PB=0
%
%   written by Try Hard
%   $Revision: 0.0.0.2 $  $Date: 2013/08/24 11:11:11 $
if nargin<1 | nargin>3 && ~isempty(cmap)
    error( 'spyc:InvalidNumArg', 'spyc takes one to three inputs')
    return
end
if isempty(sA)
    error( 'spyc:InvalidArg', 'sparse matrix is empty')
    return
end
if nargin>1 && ~isempty(cmap)
    % colorspy does not check whether your colormap is valid!
    if ~isnumeric(cmap)
        cmap=colormap(cmap);        
    end
else
    cmap=flipud(colormap('autumn'));
end
if nargin<3 || isempty(pb)
    pb=1;
end
    
indx=find(sA);
[Nx Ny]=size(sA);
sA=full(sA(indx));
ns = length(indx);
[ix iy]=ind2sub([Nx Ny],indx);
imap = round((sA-min(sA))*63/(max(sA)-min(sA)))+1;
figure, hold on
colormap(cmap)
scatter(ix,iy,[],imap,'Marker','.','SizeData',200)
set(gca,'ydir','reverse')
axis equal;
xlabel(['nz = ' num2str(ns)])
axis([0 Nx 0 Ny])
box on
if pb
    colorbar
end