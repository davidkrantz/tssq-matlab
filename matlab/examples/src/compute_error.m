function errmax = compute_error(uref1, uref2, uref3, q1, q2, q3)
    unorm = norm([uref1(:);uref2(:);uref3(:)], inf);
    err1 = abs(uref1-q1) ./ unorm;
    err2 = abs(uref2-q2) ./ unorm;
    err3 = abs(uref3-q3) ./ unorm;
    errmax = max(max(err1, err2), err3);
end
