function all_steady_state_forces()
    bulk = [1/2, 1, 2]; % the bulk moduli to be tested
    shear = [1/2, 1, 2]; % the shear moduli to be tested
    for k = bulk
        for mu = shear
            steady_state_forces(16, k, mu);
        end
    end
end