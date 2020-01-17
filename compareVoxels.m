function compareVoxels(f,struct,term,xa,xb,c)

da = f.structs{struct}.A*xa;
subplot(1,3,1)
plot(da,'.'), hold on
plot([0 length(da)],[c c])
ylabel('Dose')
title('Our Method Ax')

subplot(1,3,2)
plot(f.structs{struct}.terms{term}.w + c,'.'), hold on
plot([0 length(f.structs{struct}.terms{term}.w)],[c c])
xlabel('Index')
title('Our Method W')

db = f.structs{struct}.A*xb;
subplot(1,3,3)
plot(db,'.'), hold on
plot([0 length(db)],[c c])
title('Constraint Generation Ax')