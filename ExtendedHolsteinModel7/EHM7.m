clear; 

z = input('Enter the size of the Hamiltonian matrix'); 

A = readmatrix('h7kmap.txt'); 

h1 = zeros(z); 
h2 = zeros(z); 
h3 = zeros(z); 
h4 = zeros(z); 
h5 = zeros(z);
h6 = zeros(z); 
h7 = zeros(z); 
h8 = zeros(z); 
h9 = zeros(z); 
h10 = zeros(z);
h11 = zeros(z); 
h12 = zeros(z); 
h13 = zeros(z);
h14 = zeros(z);
h15 = zeros(z); 
h16 = zeros(z); 
h17 = zeros(z);


W = 1;                % omega -- angular frequency
t = 1;                % hopping amplitude
a = 1;                % lattice constant
i = sqrt(-1);         % imaginary i
l = 1;                % lambda
  

kValue = [0:0.1:pi]; 
k_axis = kValue/pi; 
e = zeros(numel(kValue), z);

for r = 1:numel(kValue)
    
    fj=1/((j^2+1)^(3/2));  % electron phonon interaction term
    
    k = kValue(r); 
    
    % Diagonal term

    for j= 1:z
        h1(j,j) = A(j,1)*W;
    end

    % Electron right hop -- Phonnon map left shift

    for j=1:z     
        if A(j,9) <= z           
            h2(A(j,9),j) = -t*(exp(-i*k*a));     
        end
    end
    
    % Electron left hop -- Phonnon map right shift

    for j=1:z    
        if A(j,10) <= z        
            h3(A(j,10),j) = -t*(exp(i*k*a)); 
        end
    end

    % Phonon creation

    for j=1:z    
        if A(j,11) <= z    
            h4(A(j,11),j) = -l*sqrt(A(j,2)+1);
        end
    end

    % Phonon annihilation

    for j=1:z   
        if A(j,12) <= z && A(j,12) > 0        
            h5(A(j,12),j) = -l*sqrt(A(j,2));    
        end
    end
    
    % Phonon creation at site 2 (right)
    
    for j=1:z    
        if A(j,13)<=z    
            h6(A(j,13),j) = -fj*l*sqrt(A(j,3)+1);   
        end
    end
    
    % Phonon annihilation at site 2 (left)
    
    for j=1:z
        if A(j,14)>0 && A(j,14)<=z   
            h7(A(j,14),j) = -fj*l*sqrt(A(j,3));
        end
    end
    
    % Phonon creation at site 3 (right)
    
    for j=1:z    
        if A(j,15)<=z    
            h8(A(j,15),j) = -fj*l*sqrt(A(j,4)+1);   
        end
    end
    
    % Phonon annihilation at site 3 (left)
    
    for j=1:z
        if A(j,16)>0 && A(j,16)<=z   
            h9(A(j,16),j) = -fj*l*sqrt(A(j,4));
        end
    end
    
    % Phonon creation at site 4 (right)
    
    for j=1:z    
        if A(j,17)<=z    
            h10(A(j,17),j) = -fj*l*sqrt(A(j,5)+1);   
        end
    end
    
    % Phonon annihilation at site 4 (left)
    
    for j=1:z
        if A(j,18)>0 && A(j,18)<=z   
            h11(A(j,18),j) = -fj*l*sqrt(A(j,5));
        end
    end
    
    % Phonon creation at site 5 (right)
    
    for j=1:z    
        if A(j,19)<=z    
            h12(A(j,19),j) = -fj*l*sqrt(A(j,6)+1);   
        end
    end
    
    % Phonon annihilation at site 5 (left)
    
    for j=1:z
        if A(j,20)>0 && A(j,20)<=z   
            h13(A(j,20),j) = -fj*l*sqrt(A(j,6));
        end
    end
    
    % Phonon creation at site 6 (right)
    
    for j=1:z    
        if A(j,21)<=z    
            h14(A(j,21),j) = -fj*l*sqrt(A(j,7)+1);   
        end
    end
    
    % Phonon annihilation at site 6 (left)
    
    for j=1:z
        if A(j,22)>0 && A(j,22)<=z   
            h15(A(j,22),j) = -fj*l*sqrt(A(j,7));
        end
    end
    
    % Phonon creation at site 7 (right)
    
    for j=1:z    
        if A(j,23)<=z    
            h16(A(j,23),j) = -fj*l*sqrt(A(j,8)+1);   
        end
    end
    
    % Phonon annihilation at site 7 (left)
    
    for j=1:z
        if A(j,24)>0 && A(j,24)<=z   
            h17(A(j,24),j) = -fj*l*sqrt(A(j,8));
        end
    end

    H = h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12+h13+h14+h15+h16+h17;     % total Hamiltonian matrix
    E = eig(H) ;                                                        % eigen values of the Hamiltonian
    
    e(r,:) = E; 
    
end

e= sort(real(e), 2); 
% Plot the energy eigenvalues against the normalized k-axis

for m = 1:z

    plot(k_axis, e(:,m), '-o');   % Plot each set of eigenvalues for each band
    hold on;
    
end

xlabel('k / \pi');                    
ylabel('Energy Eigenvalues');          
title('Energy Eigenvalues vs k');
grid on;
hold off;
