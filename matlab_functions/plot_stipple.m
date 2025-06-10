function [] = plot_stipple(lon,lat,stipple_field,markersize)

[X,Y] = meshgrid(lon,lat);
X = X'; Y = Y';

if nargin<4
    markersize = 1;
end

for i = 1:size(X,1)*size(X,2)
    if stipple_field(i)
        hold on; plot(X(i),Y(i),'k.','markersize',markersize)
    end
end