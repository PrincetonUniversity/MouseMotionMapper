%% Generate _by_day files

% leapout folder with all the _box_aligned_PREDICTED_2.h5 and box_data.mat

File_dir ='';
leap_dir= [File_dir 'leapout\'];
if ~isdir(leap_dir)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', leap_dir)
    return;
end

cd('leap_dir')
positions_pred_by_day = cell(size(files_by_day));
allTracks = cell(size(files_by_day));
for i = 1:size(files_by_day,1)
    for j = 1:size(files_by_day,2)
        if ~isnan(files_by_day{i,j})
            fprintf(1,['Processing dist for mouse ' num2str(i) ' day ' num2str(j) '\n']);
            nam_temp = dir(['OFT-',num2str(files_by_day{i,j},'%04.f'),'*box_aligned_PREDICTED_2.h5']);
            positions_temp = h5read(nam_temp.name, '/positions_pred');
            positions_pred_by_day{i,j}=positions_temp;

            nam_temp = dir(['OFT-',num2str(files_by_day{i,j},'%04.f'),'*.mat']);
            track_temp = load(nam_temp.name);
            allTracks{i,j}=track_temp.mouseInfo.centroidsF;
        
            fprintf(1,['No file for ' num2str(i) ' day ' num2str(j) '\n']);
            
        end
        clear nam_temp track_temp
        else
            fprintf(1,['No file for ' num2str(i) ' day ' num2str(j) '\n']);
            positions_pred_by_day{i,j}=NaN;
            allTracks{i,j}=NaN;
        end
        clear nam_temp positions_temp
    end
end

allEmbeddingValues = cell(size(files_by_day));
allOutputStats = cell(size(files_by_day));
allcGuesses = cell(size(files_by_day));
allpClusts = cell(size(files_by_day));
%%
cd()
for i = 1:size(files_by_day,1)
    for j = 1:size(files_by_day,2)
        % calculate distances from body part positions
        if ~isnan(files_by_day{i,j})
            fprintf(1,['Processing dist for mouse ' num2str(i) ' day ' num2str(j) '\n']);
            load(['Mouse_' num2str(i) '_day_' num2str(j) '_Distance_RE.mat'], 'cGuesses', 'pClusts');
            %allEmbeddingValues{i,j}=embeddingValues;
            %allOutputStats{i,j}=outputStats;
            allcGuesses{i,j}=cGuesses;
            allpClusts{i,j}=pClusts;
            clear cGuesses pClusts embeddingValues outputStats
        else
        end
    end
end

%%
cd()
for i = 1:size(files_by_day,1)
    for j = 1:size(files_by_day,2)
        % calculate distances from body part positions
        if ~isnan(files_by_day{i,j})
            fprintf(1,['Processing dist for mouse ' num2str(i) ' day ' num2str(j) '\n']);
            load(['Mouse_' num2str(i) '_day_' num2str(j) '_Distance_RE.mat'], 'embeddingValues', 'outputStats');
            allEmbeddingValues{i,j}=embeddingValues;
            allOutputStats{i,j}=outputStats;
            
            clear cGuesses pClusts embeddingValues outputStats
        else
        end
    end
end