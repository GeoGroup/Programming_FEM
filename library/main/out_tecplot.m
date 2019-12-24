function out_tecplot( fn )
% out the data ro tecplot
% Output 

fn='newtec';
ft=['tec',fn,'.dat'];
ftec=fopen(ft,'wt');
fprintf(ftec,'%s\n','TITLE     = "Fortran to Tecplot "');
fprintf(ftec,'%s\n','VARIABLES = "X"  "Y"  "Z"  "H" "PORE" '); % "HM"
fprintf(ftec,'ZONE T="GLOBAL",   N=%d,  E=%d, ZONETYPE=FEBrick\n', nn, nels);
fprintf(ftec,'DATAPACKING=BLOCK\n');
fprintf(ftec,'VARLOCATION=([1-5]=NODAL)\n');
for i=1:nn
    fprintf(ftec,'%f\n',g_coord(i,1));
end
for i=1:nn
    fprintf(ftec,'%f\n',g_coord(i,2));
end
for i=1:nn
    fprintf(ftec,'%f\n',g_coord(i,3));
end

for i=1:nn
    if disps(i)-g_coord(i,3)>0
        fprintf(ftec,'%f\n',disps(i));
    else
        fprintf(ftec,'%f\n',disps(i));
    end
end

for i=1:nn
    fprintf(ftec,'%f\n',10000*(disps(i)-g_coord(i,3)));
end

for i=1:nels
    fprintf(ftec,'%d %d %d %d %d %d %d %d\n',g_num(i,1),g_num(i,4),g_num(i,8),g_num(i,5),g_num(i,2),g_num(i,3),g_num(i,7),g_num(i,6));
end
fclose(ftec);



end

