% Parallel wrapper for the fitting algorithm
function ffit_output = parallelFit(parameter,fit_type,shaped_image,tr, submit, fit_name, ncoeffs, coeffs, tr_present,rsquared_threshold)

[dim_x, dim_y, dim_z, dim_te] = size(shaped_image);
linear_shape = reshape(shaped_image,dim_x*dim_y*dim_z,dim_te);

% voxl_pos=[];
% dlsn=0;
% for lsn=1:dim_x*dim_y*dim_z   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if sum(abs(linear_shape(lsn,:)))~=0
%         dlsn=dlsn+1;
%         rlinear_shape(dlsn,:)=linear_shape(lsn,:);
%         voxl_pos=[voxl_pos lsn];
% %         dlsn=dlsn+1;
%     end
% end
% number_voxels = dlsn;% number_voxels = dim_x*dim_y*dim_z;


abs_linear_shape=abs(linear_shape);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_linear_shape=sum(abs_linear_shape,2);
voxl_pos=find(abs_linear_shape);
linear_shape=linear_shape(voxl_pos,:);
number_voxels=length(voxl_pos);



% Preallocate for speed
if strfind(fit_type, 'user_input')
    fit_output = zeros([number_voxels 3*ncoeffs+2],'double');
elseif strfind(fit_type, 't2_exponential_plus_c')
    fit_output = zeros([number_voxels 9],'double');
else
    fit_output = zeros([number_voxels 6],'double');
end

% Break up parfor to allow for progress bar, and create progress bar
parallel_size = 200000; %10000; % 1000;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strfind(fit_type, 't1_ti') 
%     parallel_size = 10000;
% end

submit = 1; % Update not fast enough. Need progress bar to let user know that it's calculating

if submit  %Only show bar if doing the actual job
    h = waitbar(0,'Starting...','Name','Calculating fit...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
%     H = hObject;
    javaFrame = get(h,'JavaFrame');
    javaFrame.setFigureIcon(javax.swing.ImageIcon('MDACC.png'));
    
    setappdata(h,'canceling',0)
    cancel_button = false;
end

for n=1:parallel_size:number_voxels
    parfor m = n:min(n+(parallel_size-1),number_voxels)
%     for m = n:min(n+(parallel_size-1),number_voxels)   % parfor m = n:min(n+(parallel_size-1),number_voxels)
        maxNumCompThreads(2);
      workersCompThreads(m) = maxNumCompThreads;
        % iteration code here
        % note that z will be your iteration
        si = linear_shape(m,:)';
        si = cast(si,'double');
        
        fit_output(m,:) = fitParameter(parameter,fit_type,si,tr, fit_name, ncoeffs, coeffs, tr_present,rsquared_threshold);
    end
    
    
    if submit
        % check for cancel
        if getappdata(h,'canceling')
            cancel_button = true;
            break
        end
        % update waitbar each "parallel_size"th iteration
        waitbar(n/number_voxels,h,...
            sprintf('Completed %d of %d fits',n,number_voxels));
    end
end
% ffit_output = zeros(dim_x*dim_y*dim_z,dim_te);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(fit_type, 'user_input')
    ffit_output = zeros([dim_x*dim_y*dim_z 3*ncoeffs+2],'double');
elseif strfind(fit_type, 't2_exponential_plus_c')
    ffit_output = zeros([dim_x*dim_y*dim_z 9],'double');
else
    ffit_output = zeros([dim_x*dim_y*dim_z 6],'double');
end

ffit_output(voxl_pos,:)=fit_output;


if submit
    delete(h)       % DELETE the waitbar; don't try to CLOSE it.
    
    if cancel_button
        warning( 'Calculation canceled, data not valid' );
        return;
    end
end

% parfor n = 1:number_voxels
%
%     si = linear_shape(n,:)';
%     si = cast(si,'double');
%     fit_output(n,:) = fitT2(te,fit_type,si);
%
% end




