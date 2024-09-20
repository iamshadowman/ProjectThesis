clear; 

n = input('Enter the number of times to operate Hamiltonian'); 
x = 10;
a = zeros(1,x); 

kStates = {}; 
kStatesAll = {}; 

kStates{1} = -1;
%kStatesAll{1} = -1;

kStates{2} = a;
%kStatesAll{2} = a;

c1 = []; 
c2 = []; 
c3 = []; 
c4 = []; 
c5 = []; 
c6 = []; 
c7 = []; 
c8 = []; 
c9 = []; 
c10 = []; 
c11 = []; 
c12 = []; 
c13 = []; 
c14 = []; 
c15 = []; 
c16 = []; 


for j = 1: n
    
    %Electron right hop -- Phonon map left shift and check
    
    r = kStates{j+1}; 
    nr = numel(r); 
    q = zeros(1,nr); 
    
    for i = nr/2:-1:1
        if r(2*i) > 0
            if i == 1    
                q(2*i-1) = r(2*i);          
            else                    
                q(2*i-2) = r(2*i);     
            end
        end
    end
    for i = 1:nr/2
        if r(2*i-1) > 0
            q(2*i+1) = r(2*i-1);
        end
    end
    
    if length(q) ~= length(r)
        o = abs(length(q)-length(r)); 
        if mod(o,2) == 0   
            for i = 1:o
                q(end+1) = 0;
                for l = 2:length(kStates)
                    b= kStates{l}; 
                    if b ~= kStates{1}
                        b(end+1) = 0;                 
                        kStatesAll{l} = b;
                    end
                end 
                for l = 1:length(kStatesAll)
                    b = kStatesAll{l}; 
                    b(end+1) = 0;
                    kStatesAll{l} = b;
                end
            end
        elseif mod(o,2) ~= 0
            for i = 1:o  
                q(end+1) = 0;
            end
            for i = 1:o+1          
                for l = 2:length(kStates)
                    b= kStates{l}; 
                    b(end+1) = 0;
                    kStates{l} = b;
                end 
                for l = 1:length(kStatesAll)
                    b = kStatesAll{l}; 
                    if b ~= kStates{1}
                        b(end+1) = 0;                 
                        kStatesAll{l} = b;
                    end
                end
            end
        end
    end
            
    
    kStatesAll{end +1} = q; 
    
    unequal1 = true;
    for i = 1: length(kStates)
        if isequal(q, kStates{i})
            unequal1 = false; 
        end
    end
    if unequal1
        kStates{end +1} = q;
    end
    
    
    %Electron left hop -- Phonon map right shift and check
    
    t = kStates{j+1}; 
    nt = numel(t); 
    v = zeros(1,nt); 
    
    for i = nt/2:-1:1
        if t(2*i-1) > 0  
            if i == 1            
                v(2*i) = t(2*i-1); 
            else
                v(2*i-3) = t(2*i-1); 
            end
        end
    end 
    for i = 1:nt/2 
        if t(2*i) > 0
            v(2*i+2) = t(2*i);
        end
    end
    
    if length(v) ~= length(t)
        o = abs(length(v)-length(t)); 
        if mod(o,2) == 0   
            for i = 1:o
                v(end+1) = 0;
                for l = 2:length(kStates)
                    b= kStates{l}; 
                    b(end+1) = 0;
                    kStates{l} = b;
                end 
                for l = 1:length(kStatesAll)
                    b = kStatesAll{l};
                    if b ~= kStates{1}
                        b(end+1) = 0;                 
                        kStatesAll{l} = b;
                    end
                end
            end
        elseif mod(o,2) ~= 0
            for i = 1:o  
                v(end+1) = 0;
            end
            for i = 1:o+1
                for l = 2:length(kStates)
                    b= kStates{l}; 
                    b(end+1) = 0;
                    kStates{l} = b;
                end 
                for l = 1:length(kStatesAll)
                    b = kStatesAll{l}; 
                    if b ~= kStates{1}
                        b(end+1) = 0;                 
                        kStatesAll{l} = b;
                    end
                end
            end
        end
    end
    
    
    kStatesAll{end +1} = v; 
    
    unequal2 = true;
    for i = 1: length(kStates)
        if isequal(v, kStates{i})
            unequal2 = false; 
        end
    end
    if unequal2
        kStates{end +1} = v;
    end
    
    
    %Phonon creation at the site of electron  and check
    
    s=kStates{j+1}; 
    s(1) = s(1) + 1; 
    
    kStatesAll{end +1} = s; 
    
    unequal3 = true;
    for i = 1: length(kStates)
        if isequal(s, kStates{i})
            unequal3 = false; 
        end
    end
    if unequal3
        kStates{end +1} = s;
    end
    
    
    %Phonon annihilation at the site of electron and check
    
    u = kStates{j+1}; 
    u(1) = u(1) - 1; 
    
    if u(1) >= 0
        
        kStatesAll{end +1} = u;
        
        unequal4 = true;
        for i = 1:length(kStates)        
            if isequal(u, kStates{i})
                unequal4 = false;
            end   
        end
        if unequal4        
            kStates{end +1} = u;
        end
        
    elseif u(1) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 2 right to the electron and check
    
    f=kStates{j+1}; 
    f(2) = f(2) + 1; 
    
    kStatesAll{end +1} = f; 
    
    unequal5 = true;
    for i = 1: length(kStates)
        if isequal(f, kStates{i})
            unequal5 = false; 
        end
    end
    if unequal5
        kStates{end +1} = f;
    end
    
    %Phonon annihilation at the site 2 right to the electron and check
    
    g = kStates{j+1}; 
    g(2) = g(2) - 1; 
    
    if g(2) >= 0
        
        kStatesAll{end +1} = g;
        
        unequal6 = true;
        for i = 1:length(kStates)        
            if isequal(g, kStates{i})
                unequal6 = false;
            end   
        end
        if unequal6        
            kStates{end +1} = g;
        end
        
    elseif g(2) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 3 left to the electron and check
    
    h=kStates{j+1}; 
    h(3) = h(3) + 1; 
    
    kStatesAll{end +1} = h; 
    
    unequal7 = true;
    for i = 1: length(kStates)
        if isequal(h, kStates{i})
            unequal7 = false; 
        end
    end
    if unequal7
        kStates{end +1} = h;
    end
    
    %Phonon annihilation at the site 3 left to the electron and check
    
    p = kStates{j+1}; 
    p(3) = p(3) - 1; 
    
    if p(3) >= 0
        
        kStatesAll{end +1} = p;
        
        unequal8 = true;
        for i = 1:length(kStates)        
            if isequal(p, kStates{i})
                unequal8 = false;
            end   
        end
        if unequal8        
            kStates{end +1} = p;
        end
        
    elseif p(3) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 4 right to the electron and check
    
    f=kStates{j+1}; 
    f(4) = f(4) + 1; 
    
    kStatesAll{end +1} = f; 
    
    unequal9 = true;
    for i = 1: length(kStates)
        if isequal(f, kStates{i})
            unequal9 = false; 
        end
    end
    if unequal9
        kStates{end +1} = f;
    end
    
    %Phonon annihilation at the site 4 right to the electron and check
    
    g = kStates{j+1}; 
    g(4) = g(4) - 1; 
    
    if g(4) >= 0
        
        kStatesAll{end +1} = g;
        
        unequal10 = true;
        for i = 1:length(kStates)        
            if isequal(g, kStates{i})
                unequal10 = false;
            end   
        end
        if unequal10       
            kStates{end +1} = g;
        end
        
    elseif g(4) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 5 left to the electron and check
    
    h=kStates{j+1}; 
    h(5) = h(5) + 1; 
    
    kStatesAll{end +1} = h; 
    
    unequal11 = true;
    for i = 1: length(kStates)
        if isequal(h, kStates{i})
            unequal11 = false; 
        end
    end
    if unequal11
        kStates{end +1} = h;
    end
    
    %Phonon annihilation at the site 5 left to the electron and check
    
    p = kStates{j+1}; 
    p(5) = p(5) - 1; 
    
    if p(5) >= 0
        
        kStatesAll{end +1} = p;
        
        unequal12 = true;
        for i = 1:length(kStates)        
            if isequal(p, kStates{i})
                unequal12 = false;
            end   
        end
        if unequal12       
            kStates{end +1} = p;
        end
        
    elseif p(5) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 6 right to the electron and check
    
    f=kStates{j+1}; 
    f(6) = f(6) + 1; 
    
    kStatesAll{end +1} = f; 
    
    unequal13 = true;
    for i = 1: length(kStates)
        if isequal(f, kStates{i})
            unequal13 = false; 
        end
    end
    if unequal13
        kStates{end +1} = f;
    end
    
    %Phonon annihilation at the site 6 right to the electron and check
    
    g = kStates{j+1}; 
    g(6) = g(6) - 1; 
    
    if g(6) >= 0
        
        kStatesAll{end +1} = g;
        
        unequal14 = true;
        for i = 1:length(kStates)        
            if isequal(g, kStates{i})
                unequal14 = false;
            end   
        end
        if unequal14      
            kStates{end +1} = g;
        end
        
    elseif g(6) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    %Phonon creation at the site 7 left to the electron and check
    
    h=kStates{j+1}; 
    h(7) = h(7) + 1; 
    
    kStatesAll{end +1} = h; 
    
    unequal15 = true;
    for i = 1: length(kStates)
        if isequal(h, kStates{i})
            unequal15 = false; 
        end
    end
    if unequal15
        kStates{end +1} = h;
    end
    
    %Phonon annihilation at the site 7 left to the electron and check
    
    p = kStates{j+1}; 
    p(7) = p(7) - 1; 
    
    if p(7) >= 0
        
        kStatesAll{end +1} = p;
        
        unequal16 = true;
        for i = 1:length(kStates)        
            if isequal(p, kStates{i})
                unequal16 = false;
            end   
        end
        if unequal16       
            kStates{end +1} = p;
        end
        
    elseif p(7) < 0
        kStatesAll{end + 1} = kStates{1}; 
    end
    
    fprintf('Operated %d times.',j)
end

% indexing kStates

fprintf('Started indexing %d possible k states', length(kStates))

ind = -1; 
for j = 1:length(kStates)
    ind = ind + 1; 
    index(j) = ind; 
    fprintf('Indexed %d possible k states', j)
end

fprintf('Indexing completed.')

% indexing kStatesAll

fprintf('Generated %d states in total. Indexing according to the possible k states', length(kStatesAll))

for j = 1:length(kStatesAll)
    for k = 1:length(kStates)
        if kStatesAll{j} == kStates{k}
            indexAll(j) = k-1; 
        end
    end
    fprintf('Indexed %d states.', j)
end

fprintf('Indexing completed')

% calculation of number of phonons in the system

fprintf('Calculating number of phonons in a state and on electron site.')

fprintf('Calculating for %d possible k states.', length(kStates))

for j = 1:length(kStates)
    tp = 0;
    if j == 1
        tp = 0;
        p1 = 0; 
        p2 = 0; 
        p3 = 0; 
        p4 = 0;
        p5 = 0; 
        p6 = 0;
        p7 = 0; 
        
    else
        p = kStates{j}; 
        for i = 1:numel(p)
            tp = tp + p(i);   % calculate total number of phonons in the system
            p1 = p(1);        % check the number of phonons at the electron site
            p2 = p(2); 
            p3 = p(3); 
            p4 = p(4); 
            p5 = p(5);
            p6 = p(6); 
            p7 = p(7);
        
        end
    end
    tp_kStates(j) = tp;       % total number of phonons in a k state
    p1_kStates(j) = p1;       % total number of phonons at electron site of every k states
    p2_kStates(j) = p2; 
    p3_kStates(j) = p3; 
    p4_kStates(j) = p4;
    p5_kStates(j) = p5; 
    p6_kStates(j) = p6;
    p7_kStates(j) = p7; 

    fprintf('Calculated for k state %d', j)
end

fprintf('Calculation of phonons of all possible k states completed.')

fprintf('Calculating for all %d generated states.', length(kStatesAll))

for j = 1:length(kStatesAll)
    tp = 0;
    p = kStatesAll{j};
    if p == kStates{1}
        tp = 0;
        p1 = 0;
        p2 = 0; 
        p3 = 0; 
        p4 = 0; 
        p5 = 0;
        p6 = 0; 
        p7 = 0;
        
    else
        for i = 1:numel(p)
            tp = tp + p(i);      % calculate total number of phonons in the system       
            p1 = p(1);           % check the number of phonons at the electron site        
            p2 = p(2);        
            p3 = p(3);         
            p4 = p(4);         
            p5 = p(5); 
            p6 = p(6);         
            p7 = p(7); 
        end
    end
    tp_kStatesAll(j) = tp;   % total number of phonons in a k state
    p1_kStatesAll(j) = p1;   % total number of phonons at electron site of all k states
    p2_kStatesAll(j) = p2; 
    p3_kStatesAll(j) = p3;
    p4_kStatesAll(j) = p4;
    p5_kStatesAll(j) = p5;
    p6_kStatesAll(j) = p6;
    p7_kStatesAll(j) = p7;
    
    fprintf('Calculated for state %d', j)
end

fprintf('Calculation completed.')

% final columns

fprintf('Organizing final data in 6 columns and %d rows.', n)

for j = 1:n
    tp_ks(j) = tp_kStates(j+1);
    p1_ks(j) = p1_kStates(j+1); 
    p2_ks(j) = p2_kStates(j+1); 
    p3_ks(j) = p3_kStates(j+1);
    p4_ks(j) = p4_kStates(j+1); 
    p5_ks(j) = p5_kStates(j+1);
    p6_ks(j) = p4_kStates(j+1); 
    p7_ks(j) = p5_kStates(j+1);
    fprintf('Organized two columns of phonons in row %d', j)
end

e=1; 
for j = 1:length(indexAll)/16
    c1(j) = indexAll(e); 
    c2(j) = indexAll(e+1);
    c3(j) = indexAll(e+2);
    c4(j) = indexAll(e+3);
    c5(j) = indexAll(e+4); 
    c6(j) = indexAll(e+5); 
    c7(j) = indexAll(e+6);
    c8(j) = indexAll(e+7);
    c9(j) = indexAll(e+8); 
    c10(j) = indexAll(e+9); 
    c11(j) = indexAll(e+10);
    c12(j) = indexAll(e+11);
    c13(j) = indexAll(e+12);
    c14(j) = indexAll(e+13);
    c15(j) = indexAll(e+14);
    c16(j) = indexAll(e+15);
    e=e+16; 
    fprintf('Organized four columns of k states in row %d', j)
end

fprintf('Processing completed. Here is your requested output. Thank you.')

index;
indexAll;

kStates; 
kStatesAll; 

tp_kStates;
p1_kStates;

tp_kStatesAll;
p1_kStatesAll;

%table(tp_kStates', p1_kStates', index', kStates', 'VariableNames', {'total p', 'p at e site', 'k state', 'lattice pts'})
%table(tp_kStatesAll', p1_kStatesAll', indexAll', kStatesAll', 'VariableNames', {'total p', 'p at e site', 'k state', 'lattice pts'})
%table(c1', c2', c3', c4')
%table(tp_ks', p1_ks', c1', c2', c3', c4')
data=[tp_ks', p1_ks', p2_ks', p3_ks', p4_ks', p5_ks', p6_ks', p7_ks', c1', c2', c3', c4', c5', c6', c7', c8', c9', c10', c11', c12', c13', c14', c15', c16']
writematrix(data, 'h7kmap.txt')



