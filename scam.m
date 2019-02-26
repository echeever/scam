fprintf('\nStarted -- please be patient.\n\n');

fprintf('Netlist:');
type(fname)
fid = fopen(fname);
fileIn=textscan(fid,'%s %s %s %s %s %s %s');  % Read file (up to 7 items per line
[Name, N1, N2, arg3, arg4, arg5] = fileIn{:};  % Split each line into 6 columns
% Name, node1, node2, and value
fclose(fid);
nLines = length(Name);  % Number of lines in file.

N1=str2double(N1);
N2=str2double(N2);

tic
%Initialize


% Parse input file to get number of each type of element.
x=blanks(nLines);
for k1=1:nLines
    x(k1) = Name{k1}(1);  % Create string with first letter of each line
end

numZ=count(x,'R')+count(x,'L')+count(x,'C'); % Number of passive elements.
numV=count(x,'V');  % # of independent voltage sources
numO=count(x,'O');  % # of op amps
numI=count(x,'I');  % # of independent current sources
numE=count(x,'E');  % # of VCVS
numF=count(x,'F');  % # of CCCS
numG=count(x,'G');  % # of VCCS
numH=count(x,'H');  % # of CCVS
numNode = max([N1; N2]);   % Find highest node number

n=numNode;                  % number of nodes
m=numV+numO+numE+numH;      % number of voltage elements

G=cell(n,n);  [G{:}]=deal('0');    % G is nxn filled with '0'
B=cell(n,m);  [B{:}]=deal('0');
C=cell(m,n);  [C{:}]=deal('0');
D=cell(m,m);  [D{:}]=deal('0');
i=cell(n,1);  [i{:}]=deal('0');
e=cell(m,1);  [e{:}]=deal('0');
v=compose('v_%d',(1:n)');
j=cell(m,1);  [j{:}]=deal('0');

% create sub arrays (use Litovski's notation)

vsCnt = 0;      % Number of voltage sources
oaCnt = 0;      % Number of
eCnt = 0;       % Number of
hCnt = 0;       % Number of

for k1=1:nLines
    n1 = N1(k1);
    n2 = N2(k1);
    
    switch Name{k1}(1)
        % Passive element
        case {'R', 'L', 'C'} % RXXXXXXX N1 N2 VALUE
            switch Name{k1}(1)
                case 'R'
                    g=['1/' Name{k1}];
                case 'L'
                    g=['1/s/' Name{k1}];
                case 'C'
                    g=['s*' Name{k1}];
            end
            if (n1==0)
                G{n2,n2} = sprintf('%s + %s',G{n2,n2},g);
            elseif (n2==0)
                G{n1,n1} = sprintf('%s + %s',G{n1,n1},g);
            else
                G{n1,n1} = sprintf('%s + %s',G{n1,n1},g);
                G{n2,n2} = sprintf('%s + %s',G{n2,n2},g);
                G{n1,n2} = sprintf('%s - %s',G{n1,n2},g);
                G{n2,n1} = sprintf('%s - %s',G{n2,n1},g);
            end
            
        case 'V' % VXXXXXXX N1 N2 VALUE    (N1=anode, N2=cathode)
            vsCnt = vsCnt + 1;
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = '1';
            end
            if n2~=0
                C{vsCnt, n2} = '-1';
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
            end
            e{vsCnt}=Name{k1};
            j{vsCnt}=['I_' Name{k1}];
            
        case 'I' % IXXXXXXX N1 N2 VALUE  (Current N1 to N2)
            if n1~=0
                i{n1} = [i{n1} ' - ' Name{k1}];
            end
            if n2~=0
                i{n2} = [i{n2} ' + ' Name{k1}];
            end
            
        case 'O'  % % 0XXXXXXX N1 N2 N3 VALUE  (N1=+, N2=-, N3=Vout)
            n3 = str2double(arg3{k1});
            vsCnt = vsCnt + 1;
            j{vsCnt}=['I_' Name{k1}];
            if n1~=0
                C{vsCnt, n1} = [C{vsCnt, n1} ' + 1'];
            end
            if n2~=0
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            B{n3,vsCnt} = [B{n3,vsCnt} ' + 1'];
            
        case 'E'                           % VCVS
            vsCnt = vsCnt + 1;
            nc1 = str2double(arg3{k1});    % Control voltage, pos side
            nc2 = str2double(arg4{k1});    % Control voltage, neg side
             
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = [C{vsCnt, n1} ' + 1'];
            end
            if n2~=0
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            
            if nc1~=0
                C{vsCnt,nc1} = [C{vsCnt,nc1} ' - ' Name{k1}];
            end
            if nc2~=0
                C{vsCnt,nc2}= [C{vsCnt,nc2} ' + ' Name{k1}];
            end
            
            j{vsCnt}=['I_' Name{k1}];
            
            % For the CCCS we need the controlling current, which is
            % defined as the current through one of the voltage sources.  
            % Since this voltage may not have been defined yet (i.e., it 
            % comes later in the circuit definition file), we leave this
            % part of the matrix generation for later. 
            % For the CCCS there is nothing to add at this point.
        case 'F'    % CCCS FXXXXXXX N+ N- VNAM VALUE
            
        case 'G'    % VCCS GXXXXXXX N+ N- NC+ NC- VALUE
            nc1 = str2double(arg3{k1});    % Control voltage, pos side
            nc2 = str2double(arg4{k1});    % Control voltage, neg side
            g = Name{k1};
            
            % Create a string that shows if each node is ~= zero.
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
                    G{n1,nc1} =[G{n1,nc1} ' + ' g];
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
                    G{n1,nc1} =[G{n1,nc1} ' + ' g];
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
            
            
            % For the CCVS we need the controlling current which is defined as the
            % current through one of the voltage sources.  Since this voltage may not
            % have been defined yet (i.e., it comes later in the circuit definition
            % file), we leave this part of the matrix generation for later.
        case 'H'    % CCVS
            vsCnt = vsCnt + 1;
            if n1~=0
                B{n1,vsCnt} = [B{n1,vsCnt} ' + 1'];
                C{vsCnt, n1} = [C{vsCnt, n1} '+ 1'];
            end
            if n2~=0
                B{n2,vsCnt} = [B{n2,vsCnt} ' - 1'];
                C{vsCnt, n2} = [C{vsCnt, n2} ' - 1'];
            end
            j{vsCnt}=sprintf('I_%s',Name{k1});
    end
end

for k1=1:nLines
    n1 = N1(k1);
    n2 = N2(k1);
    switch Name{k1}(1)
        case 'H'
            cv = arg3{k1};  % Name of controlling voltages
            cvInd = find(contains(j,cv));  % Index of controlling voltage.
            hInd = find(contains(j,Name{k1})); % Index of CCVS (this element)
            D{hInd,cvInd}=Name{k1};
        case 'F'
            cv = arg3{k1}; % Name of controlling voltages
            cvInd = find(contains(j,cv));  % Index of controlling voltage
            if n1~=0
                B{n1,cvInd} = [B{n1,cvInd} ' + ' Name{k1}];
            end
            if n2~=0
                B{n2,cvInd} = [B{n2,cvInd} ' - ' Name{k1}];
            end
    end
end
%%

A = str2sym([G B; C D]); %Create and display x matrix
fprintf('\nThe A matrix: \n');
disp(A);

x=str2sym([v;j]);       %Create and display x matrix
fprintf('\nThe x matrix: \n');
disp(x);

z=str2sym([i;e]);       %Create and display v matrix
fprintf('\nThe z matrix:  \n');
disp(z);

syms([symvar(A), symvar(x), symvar(z)])

% The equation
fprintf('\nThe matrix equation: \n');
disp(A*x==z)

a= simplify(A\z);  % Get the solution

for i=1:length(a)  % Assign each solution to its output variable.
    eval(sprintf('%s = %s;',x(i),a(i)));
end

fprintf('\nThe solution:  \n');
disp(x==eval(x))

% Lastly assign any numeric values to symbolic variables.
for k1=1:nLines
    switch Name{k1}(1)
        % Circuit elements defined by three variable
        case {'R', 'L', 'C', 'V', 'I'}
            [num, status] = str2num(arg3{k1}); %#ok<ST2NM>
        case {'H', 'F'}
            [num, status] = str2num(arg4{k1}); %#ok<ST2NM>
        case {'E', 'G'}
            [num, status] = str2num(arg5{k1}); %#ok<ST2NM>
    end
    if status
        eval(sprintf('%s = %g;',Name{k1}, num));
    end
end

fprintf('\nElapsed time is %g seconds.\n',toc);



