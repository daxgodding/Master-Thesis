default(parisizemax, 4096M); \\ Minimum stack size needed to run all code
default(realprecision, 230); \\ Precision needed to run all code
code(N, dim, Chi1, Chi2, pol, deg, e)=
\\ Fill in 0 for Chi2 if there is only one character
{
    X=[Mod(Chi1, N), Mod(Chi2, N)]; L=[0, 0];
    \\ Creates vector with the character(s) and vector for the L-function(s)
    d=0; for(i=1, 2, if(X[i]!=0, d=d+1)); \\ Assigns d to the number of characters
    if(X[1]!=Mod(0, N), L[1]=lfunmf(M1=mfinit([N, 1, X[1]], 0), mfeigenbasis(M1)[e]));
    if(X[2]!=Mod(0, N), L[2]=lfunmf(M2=mfinit([N, 1, X[2]], 0), mfeigenbasis(M2)[e]));
    \\ Creates the space(s) of modular forms of type (1,Chi1) and (1,Chi2) of level N
    K=nfinit(pol); \\ Coefficient field
    A1=Mat([nfeltembed(K, x)~ | x <- K.zk]);
    \\ Creates the matrix with integral basis for every embedding
    if(d==1, A=A1, A=matconcat([A1; conj(A1)]));
    \\ Adds conjugate embeddings to the matrix if d is not equal to 1
    if(dim==d, Y=Mat([lfun(L[i], 0, 1) | i <- [1..d]])~);
    \\ Computes the first order derivative of all the L-functions
    if(dim!=d, Y=Mat([lfun(L[1][i], 0, 1) | i <- [1..2]])~);
    \\ Note that if d=1, and dim=2, then L[1] is a vector of two L-functions
    if(dim==4, Y=[conj(Y[2,1]); Y[1,1]; Y[2,1]; conj(Y[1,1])]);
    \\ Compatible with sequence of matrix A
    g=polredabs(algdep(exp(matsolve(A, Y)[1, 1]), deg));
    \\ Gives polynomial with reduced coefficients with exp(A^(-1)Y[1, 1]) as root
    print(g); \\ Prints the polynomial found in the previous line
    F=nfinit([nfsplitting(g), [N]]); \\ Computes the splitting field of g
    G=galoisinit(F); \\ Computes the galois group of F
    [T, n]=galoischartable(G); \\ Computes the Charactertable of G
    C=galoisconjclasses(G); \\ Computes the conjugacy classes of G
    print([[#c, permorder(c[1])] | c <- C]);
    \\ Prints order of every conjugacy class and order their elements
    for(i=1, #T,
        if(permorder(C[i][1])==2 && polsturm(galoisfixedfield(G, C[i])[1])!=0, cc=i)
    ); \\ Computes which row corresponds with the complex conjugation
    for(i=1, #T,
        if(T[1, i] == 2 && lfunconductor(lfunartin(F, G, T[, i], n))==N
            && galoischardet(G,T[, i])[cc]==-1, print([T[, i], i]))
    ); \\ Prints the row(s) of the 2-dim odd representation(s) with conductor N
};
code(23, 1, 22, 0, x, 3, 1);
code(31, 1, 30, 0, x, 3, 1);
code(39, 1, 38, 0, x, 4, 1);
code(44, 1, 21, 0, x, 3, 1);
code(47, 2, 46, 0, x^2-x-1, 5, 1);
code(52, 2, 3, 35, x^2+x+1, 9, 1);
code(55, 1, 54, 0, x, 4, 1);
code(56, 1, 13, 0, x, 4, 1);
code(57, 2, 11, 26, x^2+x+1, 9, 1);
code(59, 1, 58, 0, x, 3, 1);
code(124, 4, 67, 87, x^4-x^2+1, 24, 1);
code(133, 4, 83, 125, x^4-x^2+1, 24, 1);
code(148, 2, 105, 117, x^2+1, 24, 1);
code(229, 2, 107, 122, x^2+1, 24, 1);
code(229, 2, 107, 122, x^2+1, 24, 2);


default(parisize, 64M); \\ Beginning stack size
default(parisizemax, 800000M); \\ Highest stack size used
default(realprecision, 1500); \\ Highest precision used
N=633; \\ Level of the newform
X=[Mod(71, 633), Mod(107, 633), Mod(188, 633), Mod(266, 633)];
\\ Vector of the Dirichlet characters
L1=lfunmf(M1=mfinit([N, 1, X[1]], 0), mfeigenbasis(M1)[2]);
L2=lfunmf(M2=mfinit([N, 1, X[2]], 0), mfeigenbasis(M2)[2]);
L3=lfunmf(M3=mfinit([N, 1, X[3]], 0), mfeigenbasis(M3)[2]);
L4=lfunmf(M1=mfinit([N, 1, X[4]], 0), mfeigenbasis(M4)[2]);
\\ L-functions for spaces of type (1,Chi1),...,(1,Chi4) with level 633
L11=lfun(L1[1],0,1); \\ Derivative of L-function at 0
L12=lfun(L1[2],0,1);
L21=conj(L12); 
L22=conj(L11);
L31=lfun(L3[1],0,1);
L32=lfun(L3[2],0,1); 
L41=conj(L31); 
L42=conj(L32);
K=nfinit(y^8-y^6+y^4-y^2+1); \\ Coefficient field
A1= Mat([nfeltembed(K, w)~ | w <- K.zk]); \\ Matrix of embeddings on integral basis of K
A=matconcat([A1;conj(A1)]); \\ Adds complex conjugate of embeddings
Y=[L11, L21, L42, L31, L22, L12, L32, L41]~; \\ Compatible with sequence of A
print(g=algdep(exp(matsolve(A,Y)[1]),120));
\\ Prints polynomial with coefficients with exp(A^(-1)Y[1]) as root



default(parisize, 64M); \\ Beginning stack size
default(parisizemax, 800000M); \\ Highest stack size used
default(realprecision, 1500); \\ Highest precision used
N=1948; \\ Level of the newform
Chi=Mod(973, 1948); \\ Dirichlet character
L=lfunmf(M1=mfinit([N, 1, Chi], 0), mfeigenbasis(M1)[1]);
\\ L-functions for spaces of type (1, Chi) with level 633
L1=lfun(L[1],0,1); \\ Derivative of L-function at 0
L2=conj(L1);
L3=lfun(L[3],0,1); 
L4=conj(L3);
K=nfinit(y^4+3*y^2+1); \\ Coefficient field
A1= Mat([nfeltembed(K, w)~ | w <- K.zk]); \\ Matrix of embeddings on integral basis of K
A=matconcat([A1;conj(A1)]); \\ Adds complex conjugate embeddings
Y=[L2,L4,L1,L3]~; \\ Compatible with sequence of A
print(g=algdep(exp(matsolve(A,Y)[1]),120));
\\ Prints polynomial with coefficients with exp(A^(-1)Y[1]) as root
