% Natalia Sucher in the Kleen Lab, UCSF
% Created 1/31/2023
% Edited 2/7/2023


tic

%   1. specify directory
% EDIT THIS TO REFLECT YOUR PATH
opscea_path = '/Users/nataliasucher/Desktop/UCSF/coding/OPSCEA/';
avg_path= [opscea_path 'OPSCEADATA/avg_change_folders/'];   %path for parameters sheet


% EDIT THIS TO REFLECT THE SYMPTOMS, MODES, AND PERDUR (SECONDS BEFORE AND AFTER SYMPTOM ONSET)
% sx_input = {'lhx','rhx','lud','rud', 'lup','rup'};
% sx_input = {'lud','rud'};
sx_input = {'lex','rex'};
% sx_input = {'lhx','rhx'};
mx_input = {'2'};
perdur_input = '15';

% UNCOMMENT AND EDIT THIS IF YOU WANT TO SPECIFY PATIENT SEIZURES MANUALLY
manual_ptsz = {'EC91_03','EC96_01','EC107_01','EC133_03', 'EC166_01','EC228_03','EC229_02'}; % specify which patient and seizure
% manual_ptsz = {'EC91_03','EC133_03'}; % specify which patient and seizure
% manual_ptsz = {'EC228_03'}; % specify which patient and seizure


% LOAD RELEVANT PYTHON INFO
pe = pyenv;

% KEEP TRACK OF FOR LOOPS 
sz_count = 0;
sxmx_count = 0;



% DELETE PREVIOUS XLSX FILES
cd([avg_path '..'])

delete neurosem_pvals.xlsx
delete all_sign_change.xlsx
delete pos_sign_change.xlsx
delete neg_sign_change.xlsx
delete elec_matrix.xlsx
delete elec_w8s.xlsx

cd(opscea_path)

if str2double(pe.Version) ~= 3.9
    pyenv('Version','/Users/nataliasucher/opt/anaconda3/bin/python3.9')
end

%   1.5. for loop 
count = 0;


for sx_i = 1:length(sx_input) % for loop throughout symptoms
    sx_name = sx_input{sx_i};

    for mx_i = 1:length(mx_input) % for loop throughout modes
        mx_name = mx_input{mx_i};

        sxmx_name = [sx_name ' ' mx_name];

        sxmx_count = sxmx_count + 1;

        for ptsz_i = 1:length(manual_ptsz)

            sz_count = sz_count + 1;

            split_manual_ptsz = split(manual_ptsz{ptsz_i},'_');
            pt_name = split_manual_ptsz{1};
            sz_name = split_manual_ptsz{2};
            ptsz_name = [pt_name '-' sz_name];

            cd(opscea_path)
            [laterality, w8s_array, anat_array] = pyrunfile("activity_change.py", ["laterality", "w8s_array", "anat_array"], sxmx_input=sxmx_name, ptsz_input=ptsz_name, perdur_input=perdur_input, avg_path=avg_path, sz_count=sz_count, sxmx_count=sxmx_count, ptsz_i=ptsz_i);
            
            activity_plot(string(laterality), w8s_array, anat_array, sxmx_name, ptsz_name, perdur_input, avg_path, opscea_path, sz_count, sxmx_count, ptsz_i, pt_name, sz_name)

        end
    end 
end 

toc