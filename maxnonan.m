function m = maxnonan(x,y,d)
x(isnan(x))=-Inf;
y(isnan(y))=-Inf;
m=max(x,y,d);
end