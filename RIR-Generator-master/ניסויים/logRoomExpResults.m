function logRoomExpResults( results_table,filename,exp_disc )
% logRoomExpResults writes to a file the experiment result
writetable(results_table,'table_dummy.txt','Delimiter','\t','WriteRowNames',true)
table = fileread('table_dummy.txt');
tab_str = sprintf('\t');
table = ['   ' repmat(tab_str,1,4)  table(4:end)];
fileID = fopen(filename,'a');
fprintf(fileID,'%c',exp_disc);
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');
fprintf(fileID,'%c',table);
fprintf(fileID,'\r\n');
fprintf(fileID,'\r\n');


fclose(fileID);

end

