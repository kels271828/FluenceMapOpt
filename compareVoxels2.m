function compareVoxels2(f,struct,term,xa,xb,xc)

c = f.structs{struct}.terms{term}.dose;

da = f.structs{struct}.A*xa;
min_da = min(da);
max_da = max(da);
subplot(1,4,1)
plot(da,'.'), hold on
plot([0 length(da)],[c c])
ylabel('Dose')
title('Our Method Ax')

dw = f.structs{struct}.terms{term}.w + c;
min_dw = min(dw);
max_dw = max(dw);
subplot(1,4,2)
plot(dw,'.'), hold on
plot([0 length(f.structs{struct}.terms{term}.w)],[c c])
xlabel('Index')
title('Our Method W')

db = f.structs{struct}.A*xb;
min_db = min(db);
max_db = max(db);
subplot(1,4,3)
plot(db,'.'), hold on
plot([0 length(db)],[c c])
title('Constraint Generation Ax')

dc = f.structs{struct}.A*xc;
min_dc = min(dc);
max_dc = max(dc);
subplot(1,4,4)
plot(dc,'.'), hold on
plot([0 length(dc)],[c c])
title('conrad Ax')

min_dose = min([min_da, min_dw, min_db, min_dc]);
max_dose = max([max_da, max_dw, max_db, max_dc]);
for i = 1:3
    subplot(1,4,i)
    axis([1 length(da) min_dose max_dose])
end