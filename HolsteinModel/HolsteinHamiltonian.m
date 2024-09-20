clear; 

z = input('Enter the size of the Hamiltonian matrix'); 

A = readmatrix('kmap.txt'); 

h1 = zeros(z); 
h2 = zeros(z); 
h3 = zeros(z); 
h4 = zeros(z); 
h5 = zeros(z);

w = 1;         % omega -- angular frequency
t = 1;         % hopping amplitude
a = 1;         % lattice constant
i = sqrt(-1);  % imaginary i
l = 1;         % lambda

kValue = [0:0.1:pi]; 
k_axis = kValue/pi; 
e = zeros(numel(kValue), z);

for r = 1:numel(kValue)
    
    k = kValue(r); 
    
    % Diagonal term

    for j= 1:z
        h1(j,j) = A(j,1)*w;
    end

    % Electron right hop -- Phonnon map left shift

    for j=1:z     
        if A(j,3) <= z           
            h2(A(j,3),j) = -t*(exp(-i*k*a));     
        end
    end
    
    % Electron left hop -- Phonnon map right shift

    for j=1:z    
        if A(j,4) <= z        
            h3(A(j,4),j) = -t*(exp(i*k*a)); 
        end
    end

    % Phonon creation

    for j=1:z    
        if A(j,5) <= z    
            h4(A(j,5),j) = -l*sqrt(A(j,2)+1);
        end
    end

    % Phonon annihilation

    for j=1:z   
        if A(j,6) <= z && A(j,6) > 0        
            h5(A(j,6),j) = -l*sqrt(A(j,2));    
        end
    end

    H = h1+h2+h3+h4+h5 ;    % total Hamiltonian matrix
    E = eig(H) ;            % eigen values of the Hamiltonian
    
    h{r,1} = H; 
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
