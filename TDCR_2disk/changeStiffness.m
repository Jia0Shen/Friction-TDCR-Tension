function tdcr_ = changeStiffness(tdcr, EI)

tdcr_ = tdcr;

tdcr_.kb = EI;
nu = 0.3;
tdcr_.kt = tdcr_.kb / (1+nu);

tdcr_.Kmat = getKmat(tdcr_.kb,tdcr_.kt,tdcr_.N);
tdcr_.inv_K = getKmat(1/tdcr_.kb,1/tdcr_.kt,tdcr_.N);

end