function m = minnonan(x,y,d)
x(isnan(x))=Inf;
y(isnan(y))=Inf;
m=min(x,y,d);
end