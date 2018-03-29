

function get_basis(ord, dict, cutoff)
   
end



# # using Polynomials
# using Plots
# using DataFrames
#
# function regression(Y,Phi,lambda=0.000001)
#     #Inputs
#     #Y: vector containing exact values for the samples
#     #Phi: matrix containing values for the basis functions
#     #lambda: real - regularization parameter
#
#     #Output:
#     #vector of weights for each basis functions
#     nb_basis_fcts = size(Phi,1)
#     Q, R = qr(Phi')
#     weights = (lambda*eye(nb_basis_fcts) + R'*R) \ (Phi * Y)
#     return weights
# end
#
#
# # F(t) = sqrt(t)
# # rho(r) = exp(-r)
# #
# # function Energy_tot(x)
# #     #total energy computation for a given configuration x
# #
# #     #Input: matrix containing the configurations in row
# #     # Ntrain = size(x,1)
# #     # Nvariables = size(x,2)
# #
# #     return F.(rho.(x[:,1])+rho.(x[:,2])) +F.(rho.(x[:,2])+rho.(x[:,3])) +F.(rho.(x[:,3])+rho.(x[:,1])) #total energy
# #     # return   F.(rho.(x[:,1])+rho.(x[:,2])+rho.(x[:,3])) #site energy
# #
# # end
#
# # N = 3; #number of particles
# #
# #
# #
# # #Generate samples
# # Ntrain = 1000;
# # Ntest = 2000;
# # LeftRcuttrain = 0.5;
# # RightRcuttrain = 3;
# #
# # LeftRcuttest = 0.5;
# # RightRcuttest = 3;
# #
# #
# # r1train = (RightRcuttrain-LeftRcuttrain)*rand(Ntrain) + LeftRcuttrain; #sample variables
# # r2train = (RightRcuttrain-LeftRcuttrain)*rand(Ntrain) + LeftRcuttrain;
# # r3train = (RightRcuttrain-LeftRcuttrain)*rand(Ntrain) + LeftRcuttrain;
# #
# # r1test = (RightRcuttest-LeftRcuttest)*rand(Ntest) + LeftRcuttest; #test variables
# # r2test = (RightRcuttest-LeftRcuttest)*rand(Ntest) + LeftRcuttest;;
# # r3test = (RightRcuttest-LeftRcuttest)*rand(Ntest) + LeftRcuttest;;
# #
# # Config_train = [r1train r2train r3train]; #configuration matrix
# # Etrain = Energy_tot(Config_train) #energies for the samples
# #
# # Config_test = [r1test r2test r3test]; #configuration matrix
# # Etest_exact = Energy_tot(Config_test); #energies for the test
# #
# #
# # Deg = collect(1:10);
#
#
#
#
#
#
#
#
# a = 1; #parameter for the exponentials
#
# Poly = (["y^0","y^1","y^2","y^3","y^4","y^5","y^6","y^7","y^8","y^9","y^10","y^11"],"Polynomials");
# Rationals = (["y^0","y^(-1)","y^(-2)","y^(-3)","y^(-4)","y^(-5)","y^(-6)","y^(-7)","y^(-8)","y^(-9)","y^(-10)","y^(-11)"],"Rationals");
# Exponentials = (["y^0","exp(-$a*y)","exp(-$a*2*y)","exp(-$a*3*y)","exp(-$a*4*y)","exp(-$a*5*y)","exp(-$a*6*y)","exp(-$a*7*y)","exp(-$a*8*y)","exp(-$a*9*y)","exp(-$a*10*y)","exp(-$a*11*y)"],"Exponentials");
# Polytimesexp = (["y^0","y^1*exp(-$a*y)","y^2*exp(-$a*y)","y^3*exp(-$a*y)","y^4*exp(-$a*y)","y^5*exp(-$a*y)","y^6*exp(-$a*y)","y^7*exp(-$a*y)","y^8*exp(-$a*y)","y^9*exp(-$a*y)","y^10*exp(-$a*y)","y^11*exp(-$a*y)"],"Polynomials times exp")
# Ratiotimesexp = ["y^0","y^(-1)*exp(-$a*y)","y^(-2)*exp(-$a*y)","y^(-3)*exp(-$a*y)","y^(-4)*exp(-$a*y)","y^(-5)*exp(-$a*y)","y^(-6)*exp(-$a*y)","y^(-7)*exp(-$a*y)","y^(-8)*exp(-$a*y)","y^(-9)*exp(-$a*y)","y^(-10)*exp(-$a*y)","y^(-11)*exp(-$a*y)"],"Rationals times exp"
#
# Basis_fct_list = (Poly, Rationals, Exponentials, Polytimesexp, Ratiotimesexp);
#
# Linfty_err = zeros(length(Deg),length(Basis_fct_list))
# RMSE = zeros(length(Deg),length(Basis_fct_list))
# NB_BASIS_FCTS = zeros(length(Deg),length(Basis_fct_list));
#
#
# for (ibf,bf) in enumerate(Basis_fct_list)
#     println(Basis_fct_list[ibf][2]);
#
#     for (ideg,deg) in enumerate(Deg)
#         #generate permutation-invariant polynomials (or other functions) - uncomment the wanted basis functions
#         Perm_inv_polys = PermPolys(deg, N, Basis_fct_list[ibf][1])
#
#         #generate the basis functions
#         basis_functions = Perm_inv_polys.polys
#
#         Nb_basis_fcts = length(basis_functions)
#         NB_BASIS_FCTS[ideg,ibf] = Nb_basis_fcts
#
#         println("the maximal degree is $deg and the number of basis functions is $Nb_basis_fcts")
#
#
#         #generate the values of the basis functions in the samples points to construct
#         #the regression matrix
#         Phitrain = zeros(Nb_basis_fcts,Ntrain)
#
#         for i = 1:Ntrain
#             for j=1:Nb_basis_fcts
#                 Phitrain[j,i] = basis_functions[j](Config_train[i,:])
#             end
#         end
#
#         #compute the regression ie weights
#         Weights = regression(Etrain,Phitrain,1e-9)
#
#
#         #compute the test energy
#         Phitest = zeros(Nb_basis_fcts,Ntest)
#         for i = 1:Ntest
#             for j=1:Nb_basis_fcts
#                     Phitest[j,i] = basis_functions[j](Config_test[i,:])
#             end
#         end
#         Etest_approx = Phitest'*Weights
#
#         #compute the error
#         Linfty_err[ideg,ibf] = maximum(abs.(Etest_approx-Etest_exact))
#         RMSE[ideg,ibf] = sqrt(mean((Etest_approx-Etest_exact).*(Etest_approx-Etest_exact)))
#
#     end
# end
#

# df1 = DataFrame( deg = Deg )
# df1[Symbol("Nb_basis_fcts")] = NB_BASIS_FCTS[:,1];
# for k = 1:length(Basis_fct_list)
#     df1[Symbol(Basis_fct_list[k][2])] = RMSE[1:end,k];
# end
# println("RMSE")
#
# display(df1)


# df2 = DataFrame( deg = Deg )
# df2[Symbol("Nb_basis_fcts")] = NB_BASIS_FCTS[:,1];
# for k = 1:length(Basis_fct_list)
#     df2[Symbol(Basis_fct_list[k][2])] = Linfty_err[1:end,k];
# end
# println("Linfty_error")
#
# display(df2)
#
#
# Q = plot(Deg,RMSE[1:end,1], yscale = :log10, ylabel = "Root Mean Square Error", xlabel = "Polynomial degree",title = "Mean error versus polynomial degree",label = Basis_fct_list[1][2])
#
# for k=2:length(Basis_fct_list)
#     plot!(Deg,RMSE[1:end,k], label = Basis_fct_list[k][2])
# end
# display(Q)
#
# P = plot(NB_BASIS_FCTS[1:end,1],RMSE[1:end,1], yscale = :log10, ylabel = "Root Mean Square Error", xlabel = "Number of basis functions",title = "Mean error versus nb of basis functions",label = Basis_fct_list[1][2])
#
# for k=2:length(Basis_fct_list)
#     plot!(NB_BASIS_FCTS[1:end,1],RMSE[1:end,k], label = Basis_fct_list[k][2])
# end
# display(P)
