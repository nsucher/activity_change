% Natalia Sucher in the Kleen Lab, UCSF
% 1/31/2023


tic

%   1. specify directory
% EDIT THIS TO REFLECT YOUR PATH
opscea_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/';
avg_path= [opscea_path 'OPSCEADATA/avg_change_folders/'];   %path for parameters sheet


% EDIT THIS TO REFLECT THE SYMPTOMS, MODES, AND PERDUR (SECONDS BEFORE AND AFTER SYMPTOM ONSET)
% sx_input = {'lhx','rhx','lud','rud', 'lup','rup'};
sx_input = {'lhx','rhx'};
mx_input = {'2'};
perdur_input = '15';

% UNCOMMENT AND EDIT THIS IF YOU WANT TO SPECIFY PATIENT SEIZURES MANUALLY
% manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure
manual_ptsz = {'EC91_03','EC133_03'}; % specify which patient and seizure

% LOAD RELEVANT PYTHON INFO
pe = pyenv;

% KEEP TRACK OF FOR LOOPS 
sz_count = 0;
sxmx_count = 0;



% DELETE PREVIOUS XLSX FILES
cd(avg_path)
cd ..
if exist('Neurosemiology.xlsx')
    delete Neurosemiology.xlsx
    disp('aha! you nasty rascal.')
end

if exist('pos_sign_change.xlsx')
    delete pos_sign_change.xlsx
    disp('aha! you positive rascal.')
end

if exist('neg_sign_change.xlsx')
    delete neg_sign_change.xlsx
    disp('aha! you negative rascal.')
end


cd(opscea_path)

if str2double(pe.Version) ~= 3.9
    pyenv('Version','/Users/nataliasucher/opt/anaconda3/bin/python3.9')
end

%   1.5. for loop 
count = 0;


for sx_i = 1:length(sx_input) % for loop throughout symptoms
%     disp(count)
    sx_name = sx_input{sx_i};

    for mx_i = 1:length(mx_input) % for loop throughout modes
%         count = count + 1
        mx_name = mx_input{mx_i};


        sxmx_name = [sx_name ' ' mx_name];

        sxmx_count = sxmx_count + 1;

        if exist('manual_ptsz','var')
            for ptsz_i = 1:length(manual_ptsz)

                sz_count = sz_count + 1;

                split_manual_ptsz = split(manual_ptsz{ptsz_i},'_');
                pt_name = split_manual_ptsz{1};
                sz_name = split_manual_ptsz{2};
                ptsz_name = [pt_name '-' sz_name];

                figure; % different figure for each symptom/mode combination
                title([sxmx_name ' ' ptsz_name]);

                pyrunfile("plot_change.py", sxmx_input=sxmx_name, ptsz_input=ptsz_name, perdur_input=perdur_input, avg_path=avg_path, sz_count=sz_count, sxmx_count=sxmx_count, ptsz_i=ptsz_i);
%                 disp(toc)

%                 cd(avg_path)
% 
%                 % XLSX FILES
% 
%                 cd(opscea_path)
                
            end

        end % ~exist(manual_pt,'var')


                    

% %   2. run plot_change.py 

%         pyrunfile("plot_change.py", sxmx_input=sxmx_name, ptsz_input=ptsz_name)

%              "plot_change.py",pt='EC91',sz=);
% %       A. input 
% %           i. patient (ex.EC127)          
% %           ii. symptom (ex. rhx = right head turn)
% %           iii. mode (1 = auto, 2 = tonic, 3 = clonic)
% %           iv. perdur (single number (such as 15) that is a sec window before and after symptom onset)
% %       B. output 
% %           i. neurosemiology.xlsx 
% %   3. run sign_change.py
% %       A. input
% %           i. directory
% %       B. output 
% %           i. em_LL.xlsx
% %           ii. output sign_change.xlsx 
% %           iii. pos_sign_change.xlsx
% %           iv. neg_sign_change.xlsx
% %   4. run pval_color.m
% %       A. input
% %           i. directory
% 
% 
    end % for mx_i = 1:length(mx_input)
end % for sx_i = 1:length(sx_input)

toc