function [rect] = mostrect(mask)
	[y,x] = find(mask);
	x=round(mean(x));
	y=round(mean(y));

	[b,a] = size(mask);
	s=.1;

	rowstocolumns = a/b;

	ii=1
	rect = [x-1 y-1 2 2];
	testrect = rect;
	while imcrop(mask, testrect)==1
	    rect = testrect;
	    testrect = [x-round(ii*rowstocolumns) y-ii 2*round(ii*rowstocolumns) 2*ii];
	    ii=ii+1;
	end
	testrect
end