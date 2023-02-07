function activity_plot(laterality, w8s_array, anat_array, sxmx_name, ptsz_name, perdur_input, avg_path, opscea_path, sz_count, sxmx_count, ptsz_i, pt_name, sz_name)

cd(avg_path)
cd ..

all_T = readtable('all_sign_change.xlsx', 'Sheet', sxmx_name,'VariableNamingRule','preserve');
pos_T = readtable('pos_sign_change.xlsx', 'Sheet', sxmx_name,'VariableNamingRule','preserve');
neg_T = readtable('neg_sign_change.xlsx', 'Sheet', sxmx_name,'VariableNamingRule','preserve');
pv_T = readtable('neurosem_pvals.xlsx', 'Sheet', sxmx_name,'VariableNamingRule','preserve');
w8s_T = readtable('elec_w8s.xlsx','Sheet', sxmx_name, 'VariableNamingRule','preserve');

em_T = readtable('elec_matrix.xlsx','VariableNamingRule','preserve');


all_m = table2array(all_T(1:end,ptsz_i+1)); %extract linelength meandiff (positive or negative value of deviation from mean line length) 
pos_m = table2array(pos_T(1:end,ptsz_i+1)); %extract linelength meandiff (positive value of deviation from mean line length) 
neg_m = table2array(neg_T(1:end,ptsz_i+1)); %extract linelength meandiff (negative value of deviation from mean line length) 
pv_m = table2array(pv_T(1:end,ptsz_i+1)); % pvals per ll meandiff of neurosem across all electrodes
w8s_m = table2array(w8s_T(1:end,ptsz_i+1));
em_m = table2array(em_T(1:end,2:end));


% sign_m = {all_m,pos_m,neg_m};
sign_m = {all_m};
% sign_m = {w8s_m};



% convert py variables to matlab 

if strcmpi(laterality,sxmx_name(1)) ~= 1 %contralateral only 

    if w8s_array ~= py.NoneType
        py_w8s_cell = cell(w8s_array);
        py_anat_cell = cell(anat_array);
        
        anat_cell = {};
        w8_cell = {};
    
        for w = 1:length(py_w8s_cell)
            for c = 1:length(py_w8s_cell{w})
                anat_cell{w} = string(py_anat_cell{w}{c});
                w8_cell{w}{c} = double(py_w8s_cell{w}{c});
            end
        end
    
        for sign_i = 1:length(sign_m)
    
    
            sign_now = sign_m{1,sign_i};
            
            figure(sign_i)
            figure('Name',[sxmx_name ' ' ptsz_name]); % different figure for each symptom/mode combination
            subplot(2,1,1)
    
            cd(opscea_path)
            
            getbrain4_ns(pt_name,sz_name,1,0,laterality); %display brain 
            shading flat
            
            
            hold on;
            %SEPARATE EMNUM INTO COLUMNS OF 3 TO GET EACH ELEC MATRIX FOR EACH PTSZ
            szxyz=[];
    
            idx=(ptsz_i-1)*3;
    
            szxyz(:,:)=em_m(:,idx+1:idx+3);
    
            sz_w8s = squeeze(sign_now);
            sz_nns = ~isnan(sz_w8s);
            
            if any(pos_m)
                pos_mean_sz_LL = mean(pos_m(~isnan(pos_m(:,1))))*3.5;
            end
            if any(neg_m)
                neg_mean_sz_LL = mean(neg_m(~isnan(neg_m(:,1))))*3.5;
            end
    
            if exist("pos_mean_sz_LL","var") && exist("neg_mean_sz_LL","var")
                mean_sz_LL = mean([abs(pos_mean_sz_LL) abs(neg_mean_sz_LL)]);
            elseif exist("pos_mean_sz_LL","var") && ~exist("neg_mean_sz_LL","var")
                mean_sz_LL = pos_mean_sz_LL;
            elseif ~exist("pos_m","var") && exist("neg_mean_sz_LL","var")
                mean_sz_LL = abs(neg_mean_sz_LL);
            end
            cax = [-mean_sz_LL mean_sz_LL];
    
            if ~isnan(mean_sz_LL)
            
            %FIND ELECMATRIX FOR SZ
                ptpath = [avg_path pt_name];
                load([ptpath '/Imaging/elecs/clinical_elecs_all.mat'])
           
                %IMPORT PARAMS
                meshpath='/Imaging/Meshes/';
                Rcortex=load([ptpath meshpath pt_name '_rh_pial.mat']); 
                    loaf.rpial=Rcortex; 
                    Rcrtx=Rcortex.cortex; 
                    clear Rcortex
                Lcortex=load([ptpath meshpath pt_name '_lh_pial.mat']); 
                    loaf.lpial=Lcortex; 
                    Lcrtx=Lcortex.cortex; 
                    clear Lcortex
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                [~,prm_allPtSz]=xlsread([opscea_path 'OPSCEAparams'],'params'); 
                    fields_SZ=prm_allPtSz(1,:); % header for columns of seizure parameters
                    prm=prm_allPtSz(strcmp(pt_name,prm_allPtSz(:,1))&strcmp(sz_name,prm_allPtSz(:,2)),:);
                    if isempty(prm); error(['ATTENTION: No entry exists for ' pt_name ' seizure ' sz_name ' in the params master sheet']); end
            
                % Import parameters for patient's specific plot (layout of video frame)
                [~,plt]=xlsread([opscea_path 'OPSCEAparams'],pt_name); 
                    fields_PLOT=plt(1,:); plt(1,:)=[]; % header for columns of plotting parameters
                    plottype=plt(:,strcmpi(fields_PLOT,'plottype')); %type of plot for each subplot (accepts: iceeg, surface, depth, or colorbar)
            
                cd 
            
                surfaces=plt(:,strcmpi(fields_PLOT,'surfaces'));
                surfacesopacity=plt(:,strcmpi(fields_PLOT,'surfacesopacity'));
            
                %% plot the weights on the brain
                params_cax=str2double(regexp(prm{strcmp('cax',fields_SZ)},',','split'));         %color axis for heatmap
                params_gsp=str2double(prm{strcmp('gsp',fields_SZ)}); %gaussian spreading parameter (default 10)
                for j=1:size(plt,1) 
                          srf=regexp(surfaces{j},',','split'); % list the specific surfaces wanted for this subplot
                          srfalpha=regexp(surfacesopacity{j},',','split'); % list their corresponding opacities (values from 0 to 1; 0=invisible, 1=opaque)
                          if length(srf)~=length(srfalpha); msgbox('Number of surface to plot does not match number of alpha designations, check excel sheet'); return; end
                          acceptedterms={'rcortex','lcortex','rhipp','lhipp','ramyg','lamyg','wholebrain'};
                            for s=1:length(srf)
                                srf{s}=lower(srf{s}); %convert to lower case for easier string matching
                              if ~isempty(intersect(srf{s},acceptedterms)) %make sure user specified accepted terminologies for the meshes
                                switch srf{s} %see below for case "wholebrain"s
                                    case 'rcortex'
                                        srfplot=Rcrtx; 
                                    case 'lcortex'
                                        srfplot=Lcrtx; 
                                end
                              else
                                if laterality == 'r'
                                    srfplot=Rcrtx;
                                elseif laterality == 'l'
                                    srfplot=Lcrtx;
                                end
                              end 
                            end
        
                hh=ctmr_gauss_plot_edited(srfplot,szxyz(sz_nns,:),sz_w8s(sz_nns),cax,0,cmocean('balance'),params_gsp); 
                colorbar("southoutside",'fontsize',18)
                end
                
                lightsout
                litebrain(laterality,1)
                hold on; 
        
                for i_row=1:size(szxyz(:,1))
                    plot3(szxyz(i_row,1),szxyz(i_row,2),szxyz(i_row,3),'k.','markersize',15); 
                end 
                
    
        %%%%%%%%%%%%%%
        % % anatomy and pval plots
                subplot(2,1,2)
                hold on;
    
                for u=1:length(pv_m)   
                    if pv_m(u)<.05
                         mrkr='r*'; 
                    elseif isnan(pv_m(u))
                        continue
                    else 
                        mrkr='ko'; 
                    end
                    for w8 = 1:length(w8_cell{u})       %plot individual electrodes (w8) per neurosemiology (u)
                        plot(w8_cell{u}{w8},u*ones(length(w8_cell{u}),1),mrkr) % plot LL meandiff of each electrode
                    end
    
                end
    
                xlim(cax); 
    
                ylim([0 length(pv_m)]) 
    
                yline(0,'k-'); % vertical line at x = 0 separating positive or negative activity
                yticks(1:length(pv_m))
                yticklabels(anat_cell)
    
                for u2=1:length(pv_m)
                    if ~isnan(pv_m(u2))
                        xline(u2,'G:',.25); % horizontal lines marking neuroanatomy
                    end
                end
    
    
                alpha(1)
            end
        end
    end
end
