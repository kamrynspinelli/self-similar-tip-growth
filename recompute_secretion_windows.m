N = 128;
thresh = 0.125;
for i = 1:size(gams_l, 1)
    levels = fliplr(unique(gams_l(i,:)));
    s1_ind = find(gams_l(i,:) == max(gams_l(i,:)));
    s2_ind = find(gams_l(i,:) == max(gams_l(i,:)));
    total_gam = sum(gams_l(i,1:N/2) .* fliplr(L0Splines_l(i, end-N/2+1:end)) .* fliplr(strainls_l(i, end-N/2+1:end)));
    window_int = sum(gams_l(i,s1_ind:s2_ind-1) .* fliplr(L0Splines_l(i, end-s2_ind+2:end-s1_ind+1)) .* fliplr(strainls_l(i, end-s2_ind+2:end-s1_ind+1)));
    l = 1;
    while window_int < thresh * total_gam % assume that \gamma is either monotonic or has a single peak
        l = l + 1;
        p = find(gams_l(i,:) == levels(l));
        if p < s1_ind
            s1_ind = p;
        else
            s2_ind = p;
        end
        window_int = sum(gams_l(i,s1_ind:s2_ind-1) .* fliplr(L0Splines_l(i, end-s2_ind+2:end-s1_ind+1)) .* fliplr(strainls_l(i, end-s2_ind+2:end-s1_ind+1)));
    end
    s1 = sum(fliplr(L0Splines_l(i, end-s1_ind+2:end)) .* fliplr(strainls_l(i, end-s1_ind+2:end)));
    s2 = sum(fliplr(L0Splines_l(i, end-s2_ind+2:end)) .* fliplr(strainls_l(i, end-s2_ind+2:end)));
    W_ints_l(i,:) = [s1 s2];
    W_l(i) = s2 - s1;
end