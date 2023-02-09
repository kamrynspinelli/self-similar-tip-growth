clear all;
load('dumais-data-individual-hyphoid-fits.mat');
for prof = 1:10
    [tip_length gam_peak Z R Rtip W_int gam L0Spline strainl ksspline] = growth_anisotropic_hyphoid(fit_hyp(prof), 24);
    gam_peaks(prof) = gam_peak;
    tip_lengths(prof) = tip_length;
    ZN(prof) = Z;
    RN(prof) = R;
    R0(prof) = Rtip;
    W_ints(prof,:) = W_int;
    W(prof) = W_int(2) - W_int(1);
    gams{prof} = gam;
    L0Splines{prof} = L0Spline;
    strainls{prof} = strainl;
    kssplines{prof} = ksspline;
end
save('../cell-profiles/dumais-individual-fits-geometry-and-secretion_k_24_mu_24.mat');