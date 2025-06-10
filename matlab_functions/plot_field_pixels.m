function [] = plot_field_pixels(lon,lat,field,ctrs,cmap)

field = squeeze(field);

%colors = {[0 0 1] [0 0.5 1] [1 0.5 0] [0.75 0 0]}; %blue-red
%colors = {[0.6 0 0] [1 0.5 0] [0 0.6 1] [0 0 0.6]}; %red-blue
%colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [49 123 184]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [189 34 39]/256 [146 27 30]/256}; % alternate blue-red

%lcmap = length(ctrs)-1;
%ch = [min(ctrs) max(ctrs)];

figure;

h = pcolor(lon,lat,field'); colorbar; 
set(h,'edgecolor','none')
%caxis(ch); 
%cmap = colormap(jet(32)); 
%cmap = colormap(white0(colors,[1 1 1],32));
if nargin > 3
    caxis([ctrs(1) ctrs(end)]);
    lcmap = length(ctrs)-1;
else
    lcmap = 24;
end
if nargin > 4
    colormap(cmap);
else
    cmap = colormap(diverging0(colors,[1 1 1],lcmap));
end
if max(lat) > 60 && max(lon) > 200
    pretty_figure(600,250,'none','none',30:30:330,-90:30:90,16,{'','','90 E','','','180','','','90 W','',''},{'','60 S','30 S','EQ','30 N','60 N',''});
%elseif max(lon) > 30 && max(lon) < 40
%    pretty_figure(510,340,'none','none',-30:10:50,30:5:70,16,{'30 W','','10 W','','10 E','','30 E'},{'30 N','','40 N','','50 N','','60 N','','70 N'});
else
    pretty_figure(600,250,'none','none','none','none',16);
end