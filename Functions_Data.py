#N1, N2 funkcje kształtu 1D
shape_function1D = (lambda ksi: 0.5*(1 - ksi), \
                    lambda ksi: 0.5*(ksi + 1))
#N1, N2, N3, N4 funkcje kształtu 2D
shape_function2D = (lambda ksi, eta: 0.25*(1-ksi)*(1-eta), \
                    lambda ksi, eta: 0.25*(1+ksi)*(1-eta), \
                    lambda ksi, eta: 0.25*(1+ksi)*(1+eta), \
                    lambda ksi, eta: 0.25*(1-ksi)*(1+eta))
#dN1/dksi, dN2/dksi, dN3/dksi, dN4/dksi 2D
dN_dksi = (lambda eta: -0.25*(1-eta), \
           lambda eta: 0.25*(1-eta),  \
           lambda eta: 0.25*(1+eta),  \
           lambda eta: -0.25*(1+eta))
#dN1/deta, dN2/deta, dN3/deta, dN4/deta 2D
dN_deta = (lambda ksi: -0.25*(1-ksi), \
           lambda ksi: -0.25*(1+ksi), \
           lambda ksi: 0.25*(1+ksi),  \
           lambda ksi: 0.25*(1-ksi))
integral_points1 = [[-1./(3.**0.5), 1.], [1./(3.**0.5), 1.]]