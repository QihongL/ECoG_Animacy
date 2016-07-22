clear all; 

maxIter = 3; 
vox = 1:10; 

iter = 1; 
while iter <= maxIter
    
    
    if iter > 1 
        vox(voxUsed) = []
    end
    vox
    
    voxSel = sort(randsample(length(vox),3)')
    
    correctedVox = nan(size(voxSel));
    for v = 1 : length(voxSel)
        correctedVox(v) = voxSel(v) >= voxUsed 
    end
    
    voxUsed = voxSel
    voxUsed = sort(voxUsed)
    
    iter = iter + 1; 
end