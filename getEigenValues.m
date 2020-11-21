function [lambda1, lambda2, vect1, vect2] = getEigenValues(eigen_vects, eigen_vals)
    if eigen_vals(1, 1) > eigen_vals(2, 2)
        lambda1 = eigen_vals(1, 1);
        lambda2 = eigen_vals(2, 2);
        vect1 = eigen_vects(:, 1);
        vect2 = eigen_vects(:, 2);
    else
        lambda2 = eigen_vals(1, 1);
        lambda1 = eigen_vals(2, 2);
        vect1 = eigen_vects(:, 2);
        vect2 = eigen_vects(:, 1);
    end
end


