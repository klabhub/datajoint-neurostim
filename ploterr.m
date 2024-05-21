function varargout = ploterr(x,y,e,pv)
% Plot lines with associated error bars as shaded regions.
% x = [nrDataPoints nrLines]
% y = [nrDataPoints nrLines]
% e = [nrDataPoints nrLines] Total size of the shading (shading goes half
%       this value below y and half this value above).
% Parm/Value pairs:
% 'color'  - Matrix of colors [3 nrLines]
% 'shadingAlpha'  - Single value for all shading errors.
% lineStyle - Used for the line in the middle of the shading
% lineWidth  - Used for the line in the middle of the shading
% smooth - Smooth y with this method (see smooth) ['none']
% span - Use this span for smoothing.[5]
%
% BK - Sept 2022
arguments
    x (:,:) 
    y (:,:) 
    e (:,:) 
    pv.color = get(gca,'ColorOrder')'
    pv.shadingAlpha (1,1) double =0.8
    pv.lineStyle (1,1) string = '-';
    pv.lineWidth (1,1) double  = 1;
    pv.smooth (1,1) string = "none"
    pv.span (1,1) double = 5
end

hold on
[~,nrLines]  =size(x);
if nrLines==1 || isvector(x)
    x= repmat(x(:),[1 size(y,2)]);
end


nrY =size(y,2);
hLine = gobjects(nrY,1);
hFill = gobjects(nrY,1);
for i=1:nrY
    clrIx = mod(i-1,size(pv.color,2))+1;
    thisColor = pv.color(:,clrIx);    
    stay = ~isnan(y(:,i)) & ~isnan(e(:,i)) & ~isnan(x(:,i));
    if any(stay)
        thisX = x(stay,i);
        thisY = y(stay,i);
        thisE = e(stay,i);
        if ~strcmpi(pv.smooth,'none')
            thisY= smooth(thisX,thisY, pv.smooth,pv.span);
            thisE = smooth(thisX,thisE,pv.smooth,pv.span);
        end
        hFill(i) = fill([thisX;flipud(thisX)],[thisY-0.5*thisE;flipud(thisY+0.5*thisE)],thisColor','linestyle','none'); 
        [hFill(i).FaceAlpha] = deal(pv.shadingAlpha); 
        hold on
        hLine(i) = line(thisX,thisY,'Color',thisColor, 'LineStyle',pv.lineStyle,'LineWidth',pv.lineWidth); 
    end
end

if nargout>0
    varargout{1} = hLine;
    if nargout >1
        varargout{2} = hFill;
    end
end