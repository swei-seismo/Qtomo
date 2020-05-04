function y = glsfcn(x,transp_flag)
global G H vard varh
y = G'*((G*x)./vard) + H'*((H*x)./varh);
return