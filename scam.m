%This program takes a netlist (similar to SPICE), parses it to derive the
%circuit equations, then solves them symbolically.  
%
%Full documentation available at www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA1.html
%

if ~exist('FirstTime_rjla')
    disp(sprintf('Full documentation available at www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA1.html'));
end

disp(sprintf('\n\nStarted -- please be patient.\n'));

[Name N1 N2 arg3]=textread(fname,'%s %s %s %s ');    

tic
%Initialize
numElem=0;  %Number of passive elements.
numV=0;     %Number of independent voltage sources
numO=0;     %Number of op amps
numI=0;     %Number of independent current sources
numNode=0;  %Number of nodes, not including ground (node 0).

%Parse the input file
for i=1:length(Name),
    switch(Name{i}(1)),
        case {'R','L','C'},
            numElem=numElem+1;
            Element(numElem).Name=Name{i};
            Element(numElem).Node1=str2num(N1{i});
            Element(numElem).Node2=str2num(N2{i});
            try
                Element(numElem).Value=str2num(arg3{i});
            catch
                Element(numElem).Value=nan;
            end
        case 'V',
            numV=numV+1;
            Vsource(numV).Name=Name{i};
            Vsource(numV).Node1=str2num(N1{i});
            Vsource(numV).Node2=str2num(N2{i});
            try
                Vsource(numV).Value=str2num(arg3{i});
            catch
                Vsource(numV).Value=nan;
            end
        case 'O',
            numO=numO+1;
            Opamp(numO).Name=Name{i};
            Opamp(numO).Node1=str2num(N1{i});
            Opamp(numO).Node2=str2num(N2{i});
            Opamp(numO).Node3=str2num(arg3{i});
        case 'I'
            numI=numI+1;
            Isource(numI).Name=Name{i};
            Isource(numI).Node1=str2num(N1{i});
            Isource(numI).Node2=str2num(N2{i});
            try
                Isource(numI).Value=str2num(arg3{i});
            catch
                Isource(numI).Value=nan;
            end
    end
    numNode=max(str2num(N1{i}),max(str2num(N2{i}),numNode));
end

%Preallocate all of the cell arrays #################################
G=cell(numNode,numNode);
V=cell(numNode,1);
I=cell(numNode,1);
if ((numV+numO)~=0),
    B=cell(numNode,numV+numO);
    C=cell(numV+numO,numNode);
    D=cell(numV+numO,numV+numO);
    E=cell(numV+numO,1);
    J=cell(numV+numO,1);
end
%Done preallocating cell arrays -------------------------------------

%Fill the G matrix ##################################################
%Initially, make the G Matrix all zeros.
[G{:}]=deal('0');

%Now fill the G matrix with conductances from netlist
for i=1:numElem,
    n1=Element(i).Node1;
    n2=Element(i).Node2;
    %Make up a string with the conductance of current element.
    switch(Element(i).Name(1)),
        case 'R',
            g = ['1/' Element(i).Name];
        case 'L',
            g = ['1/s/' Element(i).Name];
        case 'C',
            g = ['s*' Element(i).Name];
    end
    
    %If neither side of the element is connected to ground
    %then subtract it from appropriate location in matrix.
    if (n1~=0) & (n2~=0),
        G{n1,n2}=[ G{n1,n2} '-' g];
        G{n2,n1}=[ G{n2,n1} '-' g];
    end
    
    %If node 1 is connected to graound, add element to diagonal
    %of matrix.
    if (n1~=0),
        G{n1,n1}=[ G{n1,n1} '+' g];
    end
    %Ditto for node 2.
    if (n2~=0),
        G{n2,n2}=[ G{n2,n2} '+' g];
    end
    
    %Go to next element.
    %     i=i+4;
end
%The G matrix is finished -------------------------------------------

%Fill the I matrix ##################################################
[I{:}]=deal('0');
for j=1:numNode,
    for i=1:numI,
        if (Isource(i).Node1==j),
            I{j}=[I{j} '-' Isource(i).Name];
        elseif (Isource(i).Node2==j),
            I{j}=[I{j} '+' Isource(i).Name];
        end
    end
end
%The I matrix is done -----------------------------------------------

%Fill the V matrix ##################################################
for i=1:numNode,
    V{i}=['v_' num2str(i)];
end
%The V matrix is finished -------------------------------------------

if ((numV+numO)~=0),
    %Fill the B matrix ##################################################
    %Initially, fill with zeros.
    [B{:}]=deal('0');
    
    %First handle the case of the independent voltage sources.
    for i=1:numV,           %Go through each independent source.
        for j=1:numNode     %Go through each node.
            if (Vsource(i).Node1==j),       %If node is first node,
                B{j,i}='1';                 %then put '1' in the matrices.
            elseif (Vsource(i).Node2==j),   %If second node, put -1.
                B{j,i}='-1';
            end
        end
    end
    
    %Now handle the case of the Op Amp
    for i=1:numO,
        for j=1:numNode
            if (Opamp(i).Node3==j),
                B{j,i+numV}='1';
            else
                B{j,i+numV}='0';
            end
        end
    end
    %The B matrix is finished -------------------------------------------
    
    
    %Fill the C matrix ##################################################
    %Initially, fill with zeros.
    [C{:}]=deal('0');
    
    %First handle the case of the independent voltage sources.
    for i=1:numV,           %Go through each independent source.
        for j=1:numNode     %Go through each node.
            if (Vsource(i).Node1==j),       %If node is first node,
                C{i,j}='1';                 %then put '1' in the matrices.
            elseif (Vsource(i).Node2==j),   %If second node, put -1.
                C{i,j}='-1';
            end
        end
    end
    
    %Now handle the case of the Op Amp
    for i=1:numO,
        for j=1:numNode
            if (Opamp(i).Node1==j),
                C{i+numV,j}='1';
            elseif (Opamp(i).Node2==j),
                C{i+numV,j}='-1';
            else
                C{i+numV,j}='0';
            end
        end
    end
    %The C matrix is finished ------------------------------------------
    
    
    %Fill the D matrix ##################################################
    %The D matrix is non-zero only for CCVS and VCVS (not included
    %in this simple implementation of SPICE)
    [D{:}]=deal('0');
    %The D matrix is finished -------------------------------------------
    
    %Fill the E matrix ##################################################
    %Start with all zeros
    [E{:}]=deal('0');
    for i=1:numV,
        E{i}=Vsource(i).Name;
    end
    %The E matrix is finished -------------------------------------------
    
    %Fill the J matrix ##################################################
    for i=1:numV,
        J{i}=['I_' Vsource(i).Name];
    end
    for i=1:numO,
        J{i+numV}=['I_' Opamp(i).Name];
    end
    %The J matrix is finished -------------------------------------------
end  %if ((numV+numO)~=0)

%Form the A, X, and Z matrices (As cell arrays of strings).
if ((numV+numO)~=0),
    Acell=[deal(G) deal(B); deal(C) deal(D)];
    Xcell=[deal(V); deal(J)];
    Zcell=[deal(I); deal(E)];
else
    Acell=[deal(G)];
    Xcell=[deal(V)];
    Zcell=[deal(I)];
end


%Declare symbolic variables #########################################
%This next section declares all variables used as symbolic variables.
%Make "s" a symbolic variable
SymString='syms s ';

%Add each of the passive elements to the list of symbolic variables.
for i=1:numElem,
    SymString=[SymString Element(i).Name ' '];
end

%Add each element of matrix J and E to the list of symbolic variables.
for i=1:numV,
    SymString=[SymString J{i} ' '];
    SymString=[SymString E{i} ' '];
end

%Add each opamp output to the list of symbolic variables.
for i=1:numO,
    SymString=[SymString J{i+numV} ' '];
end

%Add independent current sources to the list of symbolic variables.
for i=1:numI,
    SymString=[SymString Isource(i).Name ' '];
end

%Add independent voltage sources to list of symbolic variables.
for i=1:numNode,
    SymString=[SymString V{i} ' '];
end

%Evaluate the string with symbolic variables
eval(SymString);
%Done declaring symbolic variables ----------------------------------

%Create the variables A, X, and Z ###################################
%Right now the matrices Acell, Xcell and Zcell hold cell arrays of 
%strings.  These must be converted to a symbolic array.  This is
%accompplished by creating strings that represent the assignment of
%the symbolic arrays, and then evaluating these strings.

%Create assignments for three arrays
Astring='A=[';
Xstring='X=[';
Zstring='Z=[';

for i=1:length(Acell),     %for each row in the arrays.
    for j=1:length(Acell),      %for each column in matrix A.
        Astring=[Astring ' ' Acell{i,j}]; %Get element from Acell
    end
    Astring=[Astring ';'];          %Mark end of row with semicolon
    Xstring=[Xstring  Xcell{i} ';'];    %Enter element into array X;
    Zstring=[Zstring  Zcell{i} ';'];    %Enter element into array Z;
end
Astring=[Astring '];'];  %Close array assignment.
Xstring=[Xstring '];'];
Zstring=[Zstring '];'];

%Evaluate strings with array assignments.
eval([Astring ' ' Xstring ' ' Zstring])
%Done creating the variables A, X, and Z ----------------------------

%free up memory
%pack;
%Solve matrrix equation - this is the meat of the algorithm.
V=simplify(inv(A)*Z);

%Evaluate each of the unknowns in the matrix X.
for i=1:length(V),
    eval([char(X(i)) '=' char(V(i)) ';']);
end

%Assign a numeric value to each passive element, if one is provided.
for i=1:numElem
    if ~isnan(Element(i).Value),
        eval([Element(i).Name '=' num2str(Element(i).Value)  ';']);
    end
end

%Assign a numeric value to each voltage source, if one is provided.
for i=1:numV
    if ~isnan(Vsource(i).Value),
        eval([Vsource(i).Name '=' num2str(Vsource(i).Value)  ';']);
    end
end

%Assign a numeric value to each passive element, if one is provided.
for i=1:numI
    if ~isnan(Isource(i).Value),
        eval([Isource(i).Name '=' num2str(Isource(i).Value)  ';']);
    end
end

disp(sprintf('Done! Elapsed time = %g seconds.\n',toc));
beep;
disp('Netlist');
for i=1:size(Name),
    disp(sprintf('     %s %s %s %s',Name{i},N1{i},N2{i},arg3{i}));
end
disp(' ');
disp('Solved variables:');
disp(X)

if ~exist('FirstTime_rjla')
    disp(sprintf('\nFull documentation available at www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA1.html'));
    FirstTime_rjla=1;
end

