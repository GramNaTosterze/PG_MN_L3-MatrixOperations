clc
clear all
close all

% odpowiednie fragmenty kodu mozna wykonac poprzez zaznaczenie i wcisniecie F9 w Matlabie
% komentowanie/odkomentowywanie: ctrl+r / ctrl+t

% % Zadanie A
% %------------------
% N = 10;
% density = 3; % parametr decydujacy o gestosci polaczen miedzy stronami
% [Edges] = generate_network(N, density);
% %-----------------
% 
% % Zadanie B
% %------------------
% % generacja macierzy I, A, B i wektora b
% % macierze A, B i I musza byc przechowywane w formacie sparse (rzadkim)
% d = 0.85;
% b = (1 - d)/N + zeros(N,1);
% B = sparse(Edges(2,:), Edges(1,:), ones(size(Edges(1,:))), N,N );
% I = speye(N);
% L = sum(B);
% A = spdiags(1./L', 0, N,N);
% M = I - d*B*A;
% r = M\b; % Zadanie C
% %-----------------



% Zadanie D
%------------------
clc
clear all
close all

N = [500, 1000, 3000, 6000, 12000];
density = 10;


czas_Gauss = zeros(1,5);
for i = 1:5
    [Edges] = generate_network(N(i), density);
    d = 0.85;
    b = (1 - d)/N(i) + zeros(N(i),1);
    B = sparse(Edges(2,:), Edges(1,:), ones(size(Edges(1,:))), N(i),N(i) );
    I = speye(N(i));
    L = sum(B); 
    A = spdiags(1./L', 0, N(i),N(i));
    M = I - d*B*A;

    tic
    % obliczenia start
    r = M\b;
    % obliczenia stop
    czas_Gauss(i) = toc;
end
plot(N, czas_Gauss)
title("Wykres czasu obliczania układu równaniań metodą Gaussa")
xlabel("Wielkości N w kolejnych iteracjach")
ylabel("Czas obliczeń [s]")
print -dpng Wykresy/zadanieD
%------------------



% Zadanie E
%------------------
czas_Jacobi = zeros(1,5);
iter_Jacobi = zeros(1,5);
norm_Jacobi = zeros(1,1000);
norm_i = 1;
for i = 1:5
    [Edges] = generate_network(N(i), density);
    d = 0.85;
    b = (1 - d)/N(i) + zeros(N(i),1);
    B = sparse(Edges(2,:), Edges(1,:), ones(size(Edges(1,:))), N(i),N(i) );
    I = speye(N(i));
    L = sum(B); 
    A = spdiags(1./L', 0, N(i),N(i));
    M = I - d*B*A;
    
    tic
    % obliczenia start
    L = tril(M, -1);
    U = triu(M, 1);
    D = diag(diag(M));

    r = ones(N(i), 1);
    iter = 0;
    res_norm = norm(M*r - b);

    c1 = -D\(L+U);
    c2 = D\b;
    while res_norm > 10^-14
        r = c1*r+c2;
        iter = iter + 1;
        res_norm = norm(M*r - b);
        if N(i) == 1000 
            norm_Jacobi(norm_i) = res_norm;
            norm_i = norm_i + 1;
        end
    end
    % obliczenia stop
    czas_Jacobi(i) = toc;
    iter_Jacobi(i) = iter;
end

plot(N, czas_Jacobi)
title("Wykres czasu obliczania układu równaniań metodą Jacobiego")
xlabel("Wielkości N w kolejnych iteracjach")
ylabel("Czas obliczeń [s]")
print -dpng Wykresy/zadanieE_czas

plot(N, iter_Jacobi)
title("Ilość iteracji wykonanych dla danych wartośći N")
xlabel("Wielkości N w kolejnych iteracjach")
ylabel("Ilość iteracji")
print -dpng Wykresy/zadanieE_iteracje

semilogy(norm_Jacobi)
title("Wartości normalizacji dla N = 1000")
xlabel("Kolejne iteracje")
ylabel("Wartości normalizacji")
print -dpng Wykresy/zadanieE_norma
%------------------

 
% Zadanie F
%------------------
czas_GaussaSeidla = zeros(1,5);
iter_GaussaSeidla = zeros(1,5);
norm_GaussaSeidla = zeros(1,1000);
norm_i = 1;
for i = 1:5
    [Edges] = generate_network(N(i), density);
    d = 0.85;
    b = (1 - d)/N(i) + zeros(N(i),1);
    B = sparse(Edges(2,:), Edges(1,:), ones(size(Edges(1,:))), N(i),N(i) );
    I = speye(N(i));
    L = sum(B); 
    A = spdiags(1./L', 0, N(i),N(i));
    M = I - d*B*A;
    L = tril(M, -1);
    U = triu(M, 1);
    D = diag(diag(M));
    tic
    
    c1 = -(D+L);
    c2 = (D+L)\b;
    % obliczenia start
    r = ones(N(i), 1);
    iter = 0;
    res_norm = norm(M*r - b);
    while res_norm > 10^-14
        r = -c1\(U*r)+c2;
        iter = iter + 1;
        res_norm = norm(M*r - b);
        if N(i) == 1000 
            norm_GaussaSeidla(norm_i) = res_norm;
            norm_i = norm_i + 1;
        end
    end
    % obliczenia stop
    czas_GaussaSeidla(i) = toc;
    iter_GaussaSeidla(i) = iter;
end
plot(N, czas_GaussaSeidla)
title("Wykres czasu obliczania układu równaniań metodą Gaussa-Siedla")
xlabel("Wielkości N w kolejnych iteracjach")
ylabel("Czas obliczeń [s]")
print -dpng Wykresy/zadanieF_czas

title("Ilość iteracji wykonanych dla danych wartośći N")
xlabel("Wielkości N w kolejnych iteracjach")
ylabel("Ilość iteracji")
plot(N, iter_GaussaSeidla)
print -dpng Wykresy/zadanieF_iteracje

title("Wartości normalizacji dla N = 1000")
xlabel("Kolejne iteracje")
ylabel("Wartości normalizacji")
semilogy(norm_GaussaSeidla)
print -dpng Wykresy/zadanieF_norma
%------------------



% Zadanie G
%------------------
clc
clear all
close all

load('Dane_Filtr_Dielektryczny_lab3_MN.mat')
L = tril(M, -1);
U = triu(M, 1);
D = diag(diag(M));
N = size(M);
N = N(1);
% metoda Gaussa
r = M\b;
norm_Gauss = norm(M*r - b);

% metoda Jacobi
r = ones(N, 1);
res = norm(M*r - b);
a = 0;
while res > 10^-14
    r = -D\(L+U)*r+D\b;
    if res <= norm(M*r - b)
        break;
    end
    res = norm(M*r - b);
end
norm_Jacobi = res;
  % rozbieżność oznacza oddalanie się od wyniku, więc trzeba to ograniczyć
% metoda Gaussa-Seidla
r = ones(N, 1);
res = norm(M*r - b);
while res > 10^-14
    r = -(D+L)\(U*r)+(D+L)\b;
    if res <= norm(M*r - b)
        break;
    end
    res = norm(M*r - b);
end
norm_GaussSeidl = res;

print(norm_Gauss, norm_Jacobi, norm_GaussSeidl)
