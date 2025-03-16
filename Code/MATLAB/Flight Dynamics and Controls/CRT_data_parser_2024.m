%=========================================================================%
%                       CRT_data_parser.m
%
%   The m-file parses *.dat files collected using the program "crtunnel"
%   in the University of Minnesota closed return wind tunnel.  The data
%   files should be placed in a sub-directory of the directory where this
%   m-file resides.  The directory should be named "data_directory." When
%   run, this m-file will generate a matrix where each row corresponds to
%   the particular data run (angle of attack and elevator setting).  The
%   matrix has 20 columns where:
%
%   Column  1:      Row number
%   Column  2:      Angle of attack (rad)
%   Column  3:      Beta   (rad)
%   Column  4:      Elevator deflection (rad)
%   Column  5:      Rudder deflection (rad)
%   Column  6:      Aileron (rad)
%   Column  7:      Air density (kg/m^3)
%   Column  8:      Air speed (m/s)
%   Column  9:      Normal force (N)
%   Column 10:      Standard deviation of normal force (N)
%   Column 11:      Transverse force (N)
%   Column 12:      Standard deviation of transverse force (N)
%   Column 13:      Axial force (N)
%   Column 14:      Standard deviation of axial force (N)
%   Column 15:      Normal Moment (N-m)
%   Column 16:      Standard deviation of normal moment (N-m)
%   Column 17:      Transverse moment (N-m)
%   Column 18:      Standard deviation of transverse moment (N-m)
%   Column 19:      Axial moment (N-m)
%   Column 20:      Standard deviation of axial moment (N-m)
%
%   Note that columns 2 and 3 will be populated by zeros.  You must
%   populate them manually with angle of attack and elevator deflection
%   data.
%
%   Created:            1/13/2011
%   Written by:         Garrison Hoe, Demoz Gebre-Egziabher & Hamid
%                       Mokhtarzadeh
%   Last Modified:      2/18/2017 M.S. Hemati
%                       -Format 2017
%                       2/19/2018 D. Gebre-Egziabher.
%                       -Format for 2018 (Retained AOA/rudder fields
%                           in comment line).
%                       2/02/2024 Ricardo Cavalcanti Linhares.
%                       - Added beta=, elevator=, rudder=, aileron= 
%                       Also added headers.
%   Fixes
%                    
%                       2/26/2015  Peter, Sean, Demoz G. 
%                       - Format 2015
%
%                       2/11/2012 (Demoz Gebre-Egziabher)
%                       -Added GUI based file/directory selection
%                       -Increased data matrix size to 18 columns!!!!!
%                        (this will affect files that use data downstream)
%                       -Added angle of attack, elevator and rudder
%                       deflection scanning form *.dat file comment line
%
%=========================================================================%

clear all;
close all;
clc;

%   Define constants

max_elevator = 18;
d2r = pi/180;                           %   Degrees to radians
r2d = 180/pi;                           %   Radians to degrees
lb2N = 4.44822162;                      %   Pounds to Newtons
in2m = 0.0254;                          %   Inches to meters

%   Change working directory to where files are located

file_name_wild_card = '*.dat';
dialog_box_name = 'Select First File in Data Folder';
[tempFileName,file_dir_str,numFiles] = uigetfile(file_name_wild_card,dialog_box_name);

file_dir_listing_str = [file_dir_str file_name_wild_card];

file_list = dir(file_dir_listing_str);
[num_of_data_files temp] = size(file_list);

%   Sort the names so that '-2.dat' comes before '-10.dat' and so on.

temp = zeros(num_of_data_files,1);

for k=1:num_of_data_files
    
    temp1 = file_list(k).name;
    temp2 = strfind(temp1,'-');
    temp3 = strfind(temp1,'.');
    temp4 = temp1(temp2+1:temp3-1);
    temp(k) = str2num(temp4);
    
end

[new_order,old_order_index] = sort(temp,'ascend');

%   Define place holder
    
data_matrix = zeros(num_of_data_files,18);
    
%   Start cycling through each file
    
for k = 1:num_of_data_files
    
    %   Open data files for reading sequentially
    
    data_file_name=[file_dir_str file_list(old_order_index(k)).name];
    fid=fopen(data_file_name,'r');
    
    if(fid == -1)
        error(['Cannot find file named ',file_list(k).name]);
    end

    data_matrix(k,1) = k;
    
    %   Get angle of attack, elevator deflection and rudder deflection
    
    for k2 = 1:3                    %   Read lines 1 through 3
        line_read = fgets(fid);
    end
    
    user_comment_str = line_read;
    equal_sign_loc = strfind(user_comment_str,'=');

    [a1,b1] = strtok(user_comment_str(equal_sign_loc(1)+1:end),',');
%    data_matrix(k,2) = str2num(a1);         %   Beta
    data_matrix(k,3) = str2num(a1);
    

    [a1,b1] = strtok(user_comment_str(equal_sign_loc(2)+1:end),',');
    data_matrix(k,4) = str2num(a1);         %   Elevator deflection
%     
    [a1,b1] = strtok(user_comment_str(equal_sign_loc(3)+1:end),',');
    data_matrix(k,5) = str2num(a1);         %   Rudder deflection

    [a1,b1] = strtok(user_comment_str(equal_sign_loc(4)+1:end),',');
    data_matrix(k,6) = str2num(a1);         %   Aileron  
    
    
    %   Get air density

    for k2 = 4:14                   %   Read lines 4 through 11
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    data_matrix(k,7) = str2num(b2);

    %   Get air speed

    for k2 = 12:15                    %   Read lines 12 through 15
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    data_matrix(k,8) = str2num(b2);

    %   Normal Force

    for k2 = 16:22                  %   Read lines 16 through 22
        line_read = fgets(fid);
    end

    [a1,b1] = strtok(line_read,'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,9:10) = lb2N*str2num(a3);

    %   Traverse Force

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,11:12) = lb2N*str2num(a3);

    %   Axial Force

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'lb');
    data_matrix(k,13:14) = lb2N*str2num(a3);

    %   Normal Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
    data_matrix(k,15:16) = (lb2N*in2m)*str2num(a3);

    %   Traverse Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
    data_matrix(k,17:18) = (lb2N*in2m)*str2num(a3);

    %   Axial Moment

    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'in');
	data_matrix(k,19:20) = (lb2N*in2m)*str2num(a3);

    %   Angle of Attack
    [a1,b1] = strtok(fgets(fid),'=');
    [a2,b2] = strtok(b1);
    [a3,b3] = strtok(b2,'deg');
    data_matrix(k,2) = str2num(a3);
    

    fclose(fid);

end

headers = {'Run number', 'AOA', 'Beta', 'Elevator','Rudder', 'Aileron', ...
    'Density', 'Air speed', 'Normal Force','NF Std. Dev.', 'Traverse Force', 'TF Std. Dev.', 'Axial Force','AF Std. Dev.', 'Normal Moment','NM Std. Dev.', 'Traverse Moment','TM Std. Dev.', 'Axial Moment', 'Std. Dev.'};


data_matrix_with_headers = [headers; num2cell(data_matrix)];

%   Save the data file

temp = clock;
file_suffix = [num2str(temp(1)),'_', num2str(temp(2)),'_',num2str(temp(3)),...
                '_',num2str(temp(4)),'_',num2str(temp(5))];
            
save_file_name = ['CRT_data_',file_suffix,'.mat'];
eval(['save ',save_file_name,' data_matrix;']);