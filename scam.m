% This file takes as input a filename (in the variable fname) that
% describes a circuit in netlist format (similar to SPICE), and then
% performs a symbolic analysis of the circuit using Modified Nodal Analysis
% (MNA).  A full description of MNA and how to use this code is at:
% http://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA1.html
%
% Written by Erik Cheever, Swarthmore College, 2019
% echeeve1@swarthmore.edu

fprintf('\nStarted -- please be patient.\n\n');

%% Print out the netlist (a file describing the circuit with one circuit
% per line.
fprintf('Netlist:');
type(fname)
fid = fopen(fname);
fileIn=textscan(fid,'%s %s %s %s %s %s');  % Read file (up to 6 items per line
% Split each line into 6 columns, the meaning of the last 3 columns will
% vary from element to element.  The first item is always the name of the
% circuit element and the second and third items are always node numbers.
[Name, N1, N2, arg3, arg4, arg5] = fileIn{:};
% Name, node1, node2, and up to three other arguments.
fclose(fid);

nLines = length(Name);  % Number of lines in file (or elements in circuit).

N1=str2double(N1);   % Get node numbers
N2=str2double(N2);

tic                  % Begin timing.

n = max([N1; N2]);   % Find highest node number (i.e., number of nodes)

m=0; % "m" is the number of voltage sources, determined below.
for k1=1:nLines                  % Check all lines to find voltage sources
    switch Name{k1}(1)
        case {'V', 'O', 'E', 'H'}  % These are the circuit elements with
            m = m + 1;             % We have voltage source, increment m.
    end
end

% Preallocate all arrays (use Litovski's notation).
G=cell(n,n);  [G{:}]=deal('0');    % G is nxn filled with '0'
B=cell(n,m);  [B{:}]=deal('0');
C=cell(m,n);  [C{:}]=deal('0');
D=cell(m,m);  [D{:}]=deal('0');
i=cell(n,1);  [i{:}]=deal('0');
e=cell(m,1);  [e{:}]=deal('0');
j=cell(m,1);  [j{:}]=deal('0');
v=compose('v_%d',(1:n)');          % v is filled with node names

% We need to keep track of the number of voltage sources we've parsed
% so far as we go through file.  We start with zero.
vsCnt = 0;

% This loop does the bulk of filling in the arrays.  It scans line by line
% and fills in the arrays depending on the type of element found on the
% current line.
% See http://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA3.html for details.
for k1=1:nLines
    n1 = N1(k1);   % Get the two node numbers
    n2 = N2(k1);
    
    switch Name{k1}(1)
        % Passive element
        case {'R', 'L', 'C'} % RXXXXXXX N1 N2 VALUE
            switch Name{k1}(1)  % Find 1/impedance for each element type.
                case 'R'
                    g=['1/' Name{k1}];
                case 'L'
                    g=['1/s/' Name{k1}];
                case 'C'
                    g=['s*' Name{k1}];
            end
            % Here we fill in G array by adding conductance.
            % The procedure is slightly different if one of the nodes is
            % ground, so check for thos accordingly.
            if (n1==0)
                G{n2,n2} = sprintf('%s + %s',G{n2,n2},g);  % Add conductance.
            elseif (n2==0)
                G{n1,n1} = sprintf('%s + %s',G{n1,n1},g);  % Add conductance.
            else
                G{n1,n1} = sprintf('%s + %s',G{n1,n1},g);  % Add conductance.
                G{n2,n2} = sprintf('%s + %s',G{n2,n2},g);  % Add conductance.
                G{n1,n2} = sprintf('%s - %s',G{n1,n2},g);  % Sub conductance.
                G{n2,n1} = sprintf('%s - %s',G{n2,n1},g);  % Sub conductance.
            end
            
        % Independent voltage source.
        case 'V' % VXXXXXXX N1 N2 VALUE    (N1=anode, N2=cathode)
            vsCnt = vsCnt + 1;  % Keep track of which source this is.
            % Now fill in B and C arrays (again, process is slightly
            % different if one of the nodes is ground).
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = [C{vsCnt, n1} ' + 1'];
            end
            if n2~=0
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            e{vsCnt}=Name{k1};         % Add Name of source to RHS
            j{vsCnt}=['I_' Name{k1}];  % Add current through source to unknowns
            
        % Independent current source
        case 'I' % IXXXXXXX N1 N2 VALUE  (Current N1 to N2)
            % Add current to nodes (if node is not ground)
            if n1~=0
                i{n1} = [i{n1} ' - ' Name{k1}]; % subtract current from n1
            end
            if n2~=0
                i{n2} = [i{n2} ' + ' Name{k1}]; % add current to n2
            end
            
        % Op amp
        case 'O'  % 0XXXXXXX N1 N2 N3 VALUE  (N1=+, N2=-, N3=Vout)
            n3 = str2double(arg3{k1});  % This find n3
            vsCnt = vsCnt + 1;          % Keep track of number of voltage sources
            
            % Change B and C matrices as appropriate.
            B{n3,vsCnt} = [B{n3,vsCnt} ' + 1'];
            if n1~=0
                C{vsCnt, n1} = [C{vsCnt, n1} ' + 1'];
            end
            if n2~=0
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            j{vsCnt}=['I_' Name{k1}];  % Add current through source to unknowns
            
        % Voltage controlled voltage source
        case 'E'                         % VCVS
            vsCnt = vsCnt + 1;           % Keep track of number of voltage sources
            nc1 = str2double(arg3{k1});  % Control voltage, pos side
            nc2 = str2double(arg4{k1});  % Control voltage, neg side
            
            % Change B and C matrices as appropriate for output nodes.
            %  (if node is not ground)
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = [C{vsCnt, n1} ' + 1'];
            end
            if n2~=0
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            
            % Change C matrix as appropriate for input nodes
            % (if node is not ground)
            if nc1~=0
                C{vsCnt,nc1} = [C{vsCnt,nc1} ' - ' Name{k1}];
            end
            if nc2~=0
                C{vsCnt,nc2}= [C{vsCnt,nc2} ' + ' Name{k1}];
            end
            
            j{vsCnt}=['I_' Name{k1}]; % Add current through source to unknowns
            
        % Voltage controlled current source
        case 'G'    % VCCS GXXXXXXX N+ N- NC+ NC- VALUE
            nc1 = str2double(arg3{k1});    % Control voltage, pos side
            nc2 = str2double(arg4{k1});    % Control voltage, neg side
            g = Name{k1};
            
            % Create a string that shows if each node is ~= zero
            % (i.e., we find which nodes are grounded).
            myString = num2str([n1~=0, n2~=0, nc1~=0, nc2~=0]);
            myString = myString(~isspace(myString));  % Remove spaces
            % Checking all different conditions gets complicated.  There
            % may be a simpler way, but here we just brute force it and
            % check all 16 possible.
            switch myString
                case {'0000', '0011',  '0001', '0010',...
                        '0100', '1000', '1100'}
                    error('error in VCCS'); % This should never happen
                case '1111'  % All nodes are non-zero
                    G{n1,nc1} = [G{n1,nc1} ' + ' g];
                    G{n1,nc2} = [G{n1,nc2} ' - ' g];
                    G{n2,nc1} = [G{n2,nc1} ' - ' g];
                    G{n2,nc2} = [G{n2,nc2} ' + ' g];
                case '0111'  % n1 is zero - so don't include
                    G{n2,nc1} = [G{n2,nc1} ' - ' g];
                    G{n2,nc2} = [G{n2,nc2} ' + ' g];
                case '0101'
                    G{n2,nc2} = [G{n2,nc2} ' + ' g];
                case '0110'
                    G{n2,nc1} = [G{n2,nc1} ' - ' g];
                case '1011'
                    G{n1,nc1} = [G{n1,nc1} ' + ' g];
                    G{n1,nc2} = [G{n1,nc2} ' - ' g];
                case '1001'
                    G{n1,nc2} = [G{n1,nc2} ' - ' g];
                case '1010'
                    G{n1,nc1} =[G{n1,nc1} ' + ' g];
                case '1101'
                    G{n1,nc2} = [G{n1,nc2} ' - ' g];
                    G{n2,nc2} = [G{n2,nc2} ' + ' g];
                case '1110'
                    G{n1,nc1} = [G{n2,nc1} ' + ' g];
                    G{n2,nc1} = [G{n2,nc1} ' - ' g];
            end
            
        % Current controlled current source.
        % For the CCCS we need the controlling current, which is
        % defined as the current through one of the voltage sources.
        % Since this voltage may not have been defined yet (i.e., it
        % comes later in the circuit definition file), we leave this
        % part of the matrix generation for later.
        % For the CCCS there is nothing to add at this point.
        case 'F'    % CCCS FXXXXXXX N+ N- VNAM VALUE
            
        % Current controlled voltage source
        % For the CCVS we need the controlling current which is defined as the
        % current through one of the voltage sources.  Since this voltage may not
        % have been defined yet (i.e., it comes later in the circuit definition
        % file), we leave this part of the matrix generation for later.
        case 'H'    % CCVS
            vsCnt = vsCnt + 1; % Keep track of number of voltage sources
            % Change B and C as appropriate (if node is not ground)
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = [C{vsCnt, n1} '+ 1'];
            end
            if n2~=0
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            j{vsCnt}=sprintf('I_%s',Name{k1}); % Add current through source to unknowns
    end
end

% At this point all voltage sources have been parsed.  We can now go
% through and finish off the CCVS and CCCS elements (which depend on the
% current through those sources).
for k1=1:nLines
    n1 = N1(k1);
    n2 = N2(k1);
    switch Name{k1}(1)
        case 'H'
            % Here we find the indices in the matrix j:
            %    of the controlling voltage (cvInd)
            %    as well as the index of this element (hInd)
            cv = arg3{k1};  % Name of controlling voltages
            cvInd = find(contains(j,cv));  % Index of controlling voltage.
            hInd = find(contains(j,Name{k1})); % Index of CCVS (this element)
            D{hInd,cvInd}=['-' Name{k1}];  % Set the value of the D matrix.
        case 'F'
            % Here we find the index in the matrix j of the controlling
            % voltage (cvInd)
            cv = arg3{k1}; % Name of controlling voltages
            cvInd = find(contains(j,cv));  % Index of controlling voltage
            % Set the B matrix accordingly.
            if n1~=0
                B{n1,cvInd} = [B{n1,cvInd} ' + ' Name{k1}];
            end
            if n2~=0
                B{n2,cvInd} = [B{n2,cvInd} ' - ' Name{k1}];
            end
    end
end
%%  The submatrices are now complete.  Form the A, x, and z matrices,
% and solve!

A = str2sym([G B; C D]); %Create and display A matrix
fprintf('\nThe A matrix: \n');
disp(A);

x=str2sym([v;j]);       %Create and display x matrix
fprintf('\nThe x matrix: \n');
disp(x);

z=str2sym([i;e]);       %Create and display z matrix
fprintf('\nThe z matrix:  \n');
disp(z);

% Find all variables in matrices (symvar) and make them symbolic (syms)
syms([symvar(A), symvar(x), symvar(z)]);

% Displey the matrix equation
fprintf('\nThe matrix equation: \n');
disp(A*x==z)

a= simplify(A\z);  % Get the solution, this is the heart of the algorithm.

% This seems like an awkward way of doing this, if you know of a better way
% please contact me.
for i=1:length(a)  % Assign each solution to its output variable.
    eval(sprintf('%s = %s;',x(i),a(i)));
end

fprintf('\nThe solution:  \n');
disp(x==eval(x))

%% Lastly, assign any numeric values to symbolic variables.
% Go through the elements a line at a time and see if the values are valid
% numbers.  If so, assign them to the variable name.  Then you can use
% "eval" to get numberical results.
for k1=1:nLines
    switch Name{k1}(1)
        % These circuit elements defined by three variables, 2 nodes and a
        % value.  The third variable (arg3) is the value.
        case {'R', 'L', 'C', 'V', 'I'}
            [num, status] = str2num(arg3{k1}); %#ok<ST2NM>
            % Elements defined by four variables, arg4 is the value.
        case {'H', 'F'}
            [num, status] = str2num(arg4{k1}); %#ok<ST2NM>
            % Elements defined by five variables, arg5 is the value.
        case {'E', 'G'}
            [num, status] = str2num(arg5{k1}); %#ok<ST2NM>
    end
    if status  % status will be true if the argument was a valid number.
        % If the number is valid, assign it to the variable.
        eval(sprintf('%s = %g;',Name{k1}, num));
    end
end

fprintf('\nElapsed time is %g seconds.\n',toc);



