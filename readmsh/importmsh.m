function [mshfile] = importmsh(filepath)
%READMSH Generate MatLab mesh struct from .msh file from gmsh
%   Detailed explanation goes here
fid = fopen(filepath);
tline = 'start'; %placeholder to charstring to start while loop
while ischar(tline)
    % read line
    tline = fgetl(fid);
    
    % check for beginning or end of section
    if tline(1) == '$'
        section = tline(2:end);
        continue
    end
    
    % get meshformat data from .msh file
    if strcmp(section,'MeshFormat')
        mshfile.version = tline(1:3);
        if tline(5) == 1
            mshfile.file_type = "binary";
        else
            mshfile.file_type = "ASCII";
        end
        %mshfile.data_size = str2double(tline(7));
    end
    
    % get Physical names from .msh file
    if strcmp(section,'PhysicalNames')
        %how many names are there?
        num_names = str2double(tline(1));
        
        % loop over the names and write them to the struct
        for i=1:num_names
            tline = fgetl(fid); % nextline
            
            mshfile.PhysicalNames(str2double(tline(3))) = convertCharsToStrings(tline(6:end-1));
        end
    end
    
    % get entities from .msh file
    if strcmp(section,'Entities')
        nums = sscanf(tline,'%f');
        num_points = nums(1);
        num_curves = nums(2);
        num_surfaces = nums(3);
        num_volumes = nums(4);
        
       for i = 1:num_points
           tline = fgetl(fid); % nextline
           nums = sscanf(tline,'%f');
           if nums(5) > 0
               mshfile.entity_physical_tag{1}(nums(1)) = mshfile.PhysicalNames(nums(6));
           else
               mshfile.entity_physical_tag{1}(nums(1)) = "NaN";
           end
       end
       
       for i = 1:num_curves
           tline = fgetl(fid); % nextline
           nums = sscanf(tline,'%f');
           if nums(8) > 0
               mshfile.entity_physical_tag{2}(nums(1)) = mshfile.PhysicalNames(nums(9));
           else
               mshfile.entity_physical_tag{2}(nums(1)) = "NaN";
           end
       end
       
       for i = 1:num_surfaces
           tline = fgetl(fid); % nextline
           nums = sscanf(tline,'%f');
           if nums(8) > 0
               mshfile.entity_physical_tag{3}(nums(1)) = mshfile.PhysicalNames(nums(9));
           else
               mshfile.entity_physical_tag{3}(nums(1)) = "NaN";
           end
       end
       
       for i = 1:num_volumes
           tline = fgetl(fid); % nextline
           nums = sscanf(tline,'%f');
           if nums(8) > 0
               mshfile.entity_physical_tag{4}(nums(1)) = mshfile.PhysicalNames(nums(9));
           else
               mshfile.entity_physical_tag{4}(nums(1)) = "NaN";
           end
       end 
    end
    
    % get Nodes from .msh file
    if strcmp(section,'Nodes')
        nums = sscanf(tline,'%f');
        num_entity_blocks = nums(1);
        %num_nodes = nums(2);
        
        tag_index = 1;
        loc_index = 1;
        %Loop over entity blocks
        for i = 1:num_entity_blocks
            tline = fgetl(fid); % nextline
            nums = sscanf(tline,'%f');
            entity_dimension = nums(1);
            entity_index = nums(2);
            physical_tag = mshfile.entity_physical_tag{entity_dimension+1}(entity_index);
            num_nodes_in_block = nums(4);
            
            
            % write indecies and physical tags
            for j=1:num_nodes_in_block
                tline = fgetl(fid); % nextline
                num = sscanf(tline,'%f'); % readline
                mshfile.NodeIndex(tag_index) = num; % store data
                mshfile.PhysicalTag(tag_index) = physical_tag;
                tag_index = tag_index + 1;% increment index
            end
            
            % write locations
            for j=1:num_nodes_in_block
                tline = fgetl(fid); % nextline
                num = sscanf(tline,'%f')'; % readline
                mshfile.Nodes(1:2,loc_index) = num(1:2)'; % store data
                loc_index = loc_index + 1;% increment index
            end
        end
    end
    
    if strcmp(section,'Elements')
        nums = sscanf(tline,'%f');
        num_entity_blocks = nums(1);
        %num_nodes = nums(2);
        
        element_index = 1;
        bndelement_index = 1;
        %Loop over entity blocks
        for i = 1:num_entity_blocks
            tline = fgetl(fid); % nextline
            nums = sscanf(tline,'%f');
            elementtype = nums(3); %element type: 1=2-node line, 2=3-node triangle
            num_elements_in_block = nums(4);
            
            % write indecies
            for j=1:num_elements_in_block
                tline = fgetl(fid); % nextline
                num = sscanf(tline,'%f'); % readline
                if elementtype == 1 % store boundary elements
                    mshfile.BndElements(bndelement_index,:) = num(2:end);
                    bndelement_index = bndelement_index + 1;
                end
                if elementtype == 2 % store elements
                    mshfile.Elements(element_index,:) = num(2:end); % store data
                    element_index = element_index + 1;% increment index
                end
            end
        end
    end
end
fclose(fid);


% Remove redundant information
mshfile = rmfield(mshfile,'PhysicalNames');
mshfile = rmfield(mshfile,'entity_physical_tag');
mshfile = rmfield(mshfile,'NodeIndex');

mshfile.topology.element = size(mshfile.Elements,2);
mshfile.topology.boundary = size(mshfile.BndElements,2);
end

