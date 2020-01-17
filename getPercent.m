function p = getPercent(f,struct,term,x)

Ax = f.structs{struct}.A*x;
d = f.structs{struct}.terms{term}.dose;
n = f.structs{struct}.nVoxels;
p = 100*sum(Ax > d)/n;