
function val = dvhArea(f,i,x)

    dose = f.structs{i}.A*x;
    fun = @(d) 100*sum(dose > d)/f.structs{i}.nVoxels;
    val = integral(fun,0,max(dose));
    
end