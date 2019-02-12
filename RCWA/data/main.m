fileID = fopen('rect_matlab.out','w');

for i = 1 : 391
    freq = 0.05 + (i-1) * 0.005 + 0.00001;
    TRN = rcwa_hole(freq);
    fprintf(fileID, '%g\t%g\n',freq, TRN);
end

fclose(fileID);
