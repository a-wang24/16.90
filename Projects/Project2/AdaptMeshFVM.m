%This Matlab Script should be run after an initial run of an unadapted mesh
%of FVM
%This script will then calculate error indicators for each cell and the
%RefinedList vector of cell indices to be adapted
%Using the FVMrestart capability, FVM will be run again with the new
%adapted mesh from the previous solution in memory

%Compute Error Indicator
%Initiate error Indicator vector
errorInd = zeros(1,Nt);

%Calculate the difference in mach number across each interior edge and
%begin populating errorInd after multiplying with edge length
for i = 1:Ne
    iL = edge2tri(3,i);
    iR = edge2tri(4,i);
    tempErr = abs(M(iR)-M(iL))*ilengths(1,i);
    errorInd(1,iL) = errorInd(1,iL) + tempErr;
    errorInd(1,iR) = errorInd(1,iR) + tempErr;
end

%Add the contribution to error indicator from solid wall to error
%indicator vector
for i = 1:Nbe
    %index of cell adjacent to boundary edge
    iT = bedge2tri(3,i);
    %only solild wall boundaries contribute to error indicator
    if (bedge2tri(4,i) == 0)
        %find mach number normal to boundary wall
        Mn = M(i)*(dot(bnormal(:,i),[1 0]));
        tempErr = abs(Mn)*blength(1,i);
        errorInd(1,iT) = errorInd(1,iT) + tempErr;
    end
end

%Sort the error indicator vector in descending order using the sort
%function which returns the vector sorted and their original indices
%within the errorInd vector which is useful because that is their cell
%index
%Note: errorInd(OrigInd) is equal to sortedError
[sortedError,origInd] = sort(errorInd, 'descend');

%OrigInd is a vector of cell indexes sorted so the first index has the
%highest error indicator value
%Choose to use eta=.25 and adapt the meshes of the first quarter of the
%cell indexes within OrigInd
%Assign the first 25% of indexes to RefineList to use the script
%cyl_adaptmesh
eta = 1/4;
remainder = rem(length(origInd),eta);
lastInd = (length(origInd)*eta)+(remainder*eta); %ensures index is int
RefineList = origInd(1:lastInd);
RefineList = RefineList';

cyl_adaptmesh

%use restart capability of FVM
%increase ntol
ntol = n+100;
FVMrestart = 1;
FVM



