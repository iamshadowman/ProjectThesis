clear; 

z = input('Enter the size of the Hamiltonian matrix'); 

A = readmatrix('h3kmap.txt'); 

h1 = zeros(z); 
h2 = zeros(z); 
h3 = zeros(z); 
h4 = zeros(z); 
h5 = zeros(z);
h6 = zeros(z); 
h7 = zeros(z); 
h8 = zeros(z); 
h9 = zeros(z); 

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
        if A(j,5) <= z           
            h2(A(j,5),j) = -t*(exp(-i*k*a));     
        end
    end
    
    % Electron left hop -- Phonnon map right shift

    for j=1:z    
        if A(j,6) <= z        
            h3(A(j,6),j) = -t*(exp(i*k*a)); 
        end
    end

    % Phonon creation

    for j=1:z    
        if A(j,7) <= z    
            h4(A(j,7),j) = -l*sqrt(A(j,2)+1);
        end
    end

    % Phonon annihilation

    for j=1:z   
        if A(j,8) <= z && A(j,8) > 0        
            h5(A(j,8),j) = -l*sqrt(A(j,2));    
        end
    end
    
    % Phonon creation at site 2 (right)
    
    for j=1:z    
        if A(j,9)<=z    
            h6(A(j,9),j) = -fj*l*sqrt(A(j,3)+1);   
        end
    end
    
    % Phonon annihilation at site 2 (left)
    
    for j=1:z
        if A(j,10)>0 && A(j,10)<=z   
            h7(A(j,10),j) = -fj*l*sqrt(A(j,3));
        end
    end
    
    % Phonon creation at site 3 (right)
    
    for j=1:z    
        if A(j,11)<=z    
            h8(A(j,11),j) = -fj*l*sqrt(A(j,4)+1);   
        end
    end
    
    % Phonon annihilation at site 3 (left)
    
    for j=1:z
        if A(j,12)>0 && A(j,12)<=z   
            h9(A(j,12),j) = -fj*l*sqrt(A(j,4));
        end
    end

    H = h1+h2+h3+h4+h5+h6+h7+h8+h9 ;    % total Hamiltonian matrix
    E = eig(H) ;                        % eigen values of the Hamiltonian
    
    e(r,:) = E; 
    
end

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
